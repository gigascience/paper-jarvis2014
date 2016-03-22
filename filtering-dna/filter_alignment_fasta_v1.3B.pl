#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

my ($input,$windows,$min_identity,$len_specific,$help);
GetOptions(
	"input:s"=>\$input,
	"window:i"=>\$windows,
	"minidentity:f"=>\$min_identity,
	"len:i"=>\$len_specific,
	"help"=>\$help,
);
if(!$input or $help){usage();exit;}
if(!open IN1,"$input"){print "Error:The input file can not be open\n";usage();}
my $base_name=basename($input);
$base_name=~s/\.\w+$//;
$windows||=15;
$min_identity||=0.3;
$len_specific||=6;
my %hash_fa;
my %hash_identity;
my $id='';
####read fa file
while(<IN1>){
	if(/>(\S+)/){$id=$1;next;}
	chomp;
	s/\s//g;
	tr/[atcgn]/[ATCGN]/;
	$hash_fa{$id}.=$_;
}
close IN1;
my %hash_out1;
####### get the specific block
my @id=sort keys %hash_fa;
for(my $i=0; $i<@id; $i++){
	my $line='';
	for(my $k=0;$k<length($hash_fa{$id[$i]});$k++){
		my $l_tmp=substr($hash_fa{$id[$i]},$k,1);
		if($l_tmp eq 'N'){if(special_Ns($id[$i],$k,1)==1){$l_tmp=~tr/N/n/;}}
		else{
			if(calculate_identity($id[$i],$k,1)<$min_identity){$l_tmp=~tr/[ATCG]/[atcg]/;}
		}
		$line.=$l_tmp;
	}
	while($line=~/([atcg]{$len_specific,})/g){my $L=length($1);my $a="n"x$L;$line=~s/$1/$a/;}
	$line=~tr/[atcg]/[ATCG]/;
	$hash_out1{$id[$i]}=$line;
}

######### get the low identity block. The window size is $windows. If meet '---' or 'n' or 'N', the window will be forward terminated.
my %hash_out2;
for(my $i=0;$i<@id;$i++){
	my $line=$hash_out1{$id[$i]};
	my $line_tmp='';
	my @pos;$pos[0]=$pos[1]=$pos[2]=$pos[3]=0;##$pos[0] the start pos; $pos[1] the length; $pos[2] the flag 0:the prior is end by un'n', 1: is n.
	
	for(my $k=0;$k<length($line);$k++){
		my $base=substr($line,$k,1);
		$line_tmp.=$base;
		$pos[1]++;
		$pos[0]=$k if($pos[0]==0&& $pos[3]>0);
		if($base=~/n|N/){
			if($pos[2]==1 && $pos[1]<4){
				$line_tmp=~tr/[ATCG]/[atcg]/;
			}
			if(calculate_identity($id[$i],$pos[0],$pos[1])<$min_identity){
				$line_tmp=~tr/[ATCG]/[atcg]/;
			}
			$hash_out2{$id[$i]}.=$line_tmp;
			$pos[0]=$pos[1]=0;$line_tmp='';$pos[2]=1;$pos[3]++;
			next;
		}
		
		if($line_tmp=~/-{6,}/ or $k+1==length($line)){
			if(calculate_identity($id[$i],$pos[0],$pos[1])<$min_identity){$line_tmp=~tr/[ATCG]/[atcg]/;}
			$hash_out2{$id[$i]}.=$line_tmp;
			$pos[0]=$pos[1]=$pos[2]=0;$line_tmp='';$pos[3]++;
			next;
		}
		if(($pos[1]-$windows)==int(0.5*$windows+0.5)){
			$k=$k-int(0.5*$windows+0.5);$pos[3]++;
			$pos[1]=$pos[1]-int(0.5*$windows+0.5);
			my $flag_pos=$pos[1];
			#print "$k\t$pos[0]\t$pos[1]\n";###for test
			my $li;
			my $identity_0=calculate_identity($id[$i],$pos[0],$pos[1]);
			my $identity_tmp=$identity_0;
			if($identity_0<$min_identity){

				#while(1){
				#	$pos[1]--;$k--;
				#	my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1]);
				#	if($identity_1>$identity_0){$identity_0=$identity_1;}
				#	else{last;}
				#}
				#$pos[1]++;$k++;
				my $pos_tmp=$pos[1];
				if($pos_tmp==$flag_pos){
					$identity_0=$identity_tmp;
					while(1){
						$pos[1]++;$k++;
						my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1]);
						if($identity_1<$identity_0){$identity_0=$identity_1;}
						else{last;}
					}
					$pos[1]--;$k--;
				}
				$li=substr($line,$pos[0],$pos[1]);
				$li=~tr/[ATCG]/[atcg]/;
			}
			else{
				while(1){
					$pos[1]--;$k--;
					my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1]);
					if($identity_1<$identity_0){$identity_0=$identity_1;}
					else{last;}
				}
				$pos[1]++;$k++;
				my $pos_tmp=$pos[1];
				if($pos_tmp==$flag_pos){
					$identity_0=$identity_tmp;
					while(1){
						$pos[1]++;$k++;
						my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1]);
						if($identity_1>$identity_0){$identity_0=$identity_1;}
						else{last;}
					}
					$pos[1]--;$k--;
				}
				$li=substr($line,$pos[0],$pos[1]);
			}
			#print "$k\t$pos[0]\t$pos[1]\n";###test
			$hash_out2{$id[$i]}.=$li;
			$pos[0]=$pos[1]=$pos[2]=0;$line_tmp='';
			next;
		}
		
	}
	if($line_tmp){
		$pos[3]++;
		if(calculate_identity($id[$i],$pos[0],$pos[1])<$min_identity){$line_tmp=~tr/[ATCG]/[atcg]/;}
		$hash_out2{$id[$i]}.=$line_tmp;
	}
}

############# put out 
print_pos(\%hash_out2,"$base_name-$windows-$min_identity.stat.xls");
put_out_fa(\%hash_out2,"$base_name-$windows-$min_identity.filter");

############### sub 

###compared with mini_identity
sub special_Ns{
	my ($id,$start,$length)=@_;
	my $f=1;
	my $l1=substr($hash_fa{$id},$start,$length);
	if($l1 ne 'N'){$f=0;}
	else{
		foreach(@id){
	                if($_ eq $id){next;}
        	        my $l2=substr($hash_fa{$_},$start,$length);
	                if($l2 ne '-'){$f=0;}
	        }
	}
	return $f;
}

sub calculate_identity{
	my ($id,$start,$length)=@_;
	my $f=0;
	if($length==0){return $f;}
	my $l1=substr($hash_fa{$id},$start,$length);
	foreach(@id){
		if($_ eq $id){next;}
		my $l2=substr($hash_fa{$_},$start,$length);
		my $iden=get_identity($l1,$l2);
		if($iden>$f){$f=$iden;}
	}
	return $f;
}

#### put out xls file
sub print_pos{
	my ($hash,$out)=@_;
	my %tmp_hash;
	open OUTPOS,">$out";
	print OUTPOS "##the position is 1-based\n##the low identity means identity<$min_identity,the windows size is $windows\nspecies_name\tstart(gap)\tend(gap)\tstart(nogap)\tend(nogap)\tlength\ttype\n";
	foreach(sort keys %{$hash}){
		my $line=$hash->{$_};
		my $id_key=$_;
		$line=~tr/[atcg]/+/;
		my ($start1,$start2,$s1,$s2,$flag1,$flag2,$n)=(0,0,0,0,0,0,0);
		for(my $k=0;$k<length($line);$k++){
			my $a=substr($line,$k,1);
			if($a ne '-'){$n++;}
			if($flag1==0 && $a eq '+'){$flag1=1;$start1=$k+1;$start2=$n+1;}
			if($flag1==1 && $a ne '+'){
				my $l=$k-$start1+1;
				my $end2=$start2+$l-1;
				$flag1=0;
				print OUTPOS "$id_key\t$start1\t$k\t$start2\t$end2\t$l\tlow identity\n";
				$start1=$start2=0;
			}
			if($flag2==0 && $a eq 'n'){$flag2=1;$s1=$k+1;$s2=$n+1;}
			if($flag2==1 && $a ne 'n'){
				my $l=$k-$s1+1;
				my $end2=$start2+$l-1;
				$flag2=0;
				print OUTPOS "$id_key\t$s1\t$k\t$s2\t$end2\t$l\tspecific\n";
				$s1=$s2=0;
			}
		}
	}
}

###put out the changed fa file. specific block masked with 'n', low identity masked wtih lower case forms
sub put_out_fa{
	my ($hash,$out)=@_;
	if($out){open OUT,">$out";}
	foreach my $id(sort keys %{$hash}){
		if($out){print OUT ">$id\n$hash->{$id}\n";}
		else{print ">$id\n$hash->{$id}\n";}
		}
	if($out){close OUT;}
}

sub get_identity{
	my ($line1,$line2)=@_;
	my $l1=length($line1);
	my $l2=length($line2);
	my $equal_num=0;
	my $gap_num=0;
	die "the length is unequal\n" if($l1!=$l2);
	for(my $i=0;$i<$l1;$i++){
		my $a=substr($line1,$i,1);
		my $b=substr($line2,$i,1);
		if($a eq "-"){$gap_num++;next;}
		elsif($a eq $b){ $equal_num++;}
	}
	my $result=($l1==$gap_num)?1:$equal_num/($l1-$gap_num);
	return $result;
}

sub usage{
	print <<USAGE;
	This program is used to mask the specific and low identity sequences of multiple sequence alignment in FASTA format by 'n' and lowercase respectively.
	
	perl $0 --input <input.fa> --window <int number> --minidentity <float number> --len <int number>;
	or perl $0 -i <input.fa> -w <int number> -m <float number> -l <int number>;

	--input(i): multiple sequence alignment in FASTA format;
	--window(w): int number,default 15;
	--minidentity(m):  minimum identity, default 0.3;
	--len(l): the minimum length of specific, default 6;
	--help(h): print the help information on screen.
	
	The output file: *.fliter:alignment file after masking and *.xls: filtering information;
	
USAGE
}
