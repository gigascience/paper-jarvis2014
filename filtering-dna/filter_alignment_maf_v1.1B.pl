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
if(!open IN,"$input"){print "Error:The input file can not be open\n";usage();}
my $base_name=basename($input);
$base_name=~s/\.\w+$//;
open OUT,">$base_name.filter";
$windows||=15;
$min_identity||=0.3;
$len_specific||=6;

my %hash_identity;
my $id='';
my $block=0;
####read maf file#####
$/="a score=";
chomp(my $title=<IN>);
print OUT $title;
while(<IN>){
	my %hash_inf;
	my %hash_maf;
	my %hash_out;
	my %hash_mask;
	chomp;
	my @lines=split(/\n/,$_);
	foreach my $line (@lines){
		if($line=~/^s/){
			my $l=$line;
			my @inf=split(/\s+/,$line);
			$inf[6]=~tr/[atcgn]/[ATCGN]/;
			$hash_maf{$inf[1]}=$inf[6];
			$l=~s/\S+$//;
			$hash_inf{$inf[1]}=$l;
		#	print "$hash_inf{$inf[1]}\n$hash_maf{$inf[1]}\n";
		}
	
	}
	specific_block(\%hash_maf,\%hash_out);
	low_identity(\%hash_out,\%hash_mask);
	print OUT "a score=$lines[0]\n";
	foreach my $x (sort keys %hash_maf){
		print OUT "$hash_inf{$x}$hash_mask{$x}\n";
	}
	print  OUT "\n";
	undef %hash_inf;
	undef %hash_maf;
	undef %hash_out;
	undef %hash_mask;
}

####### get the specific block
sub specific_block{
	my ($hash_fa,$hash_out1)=@_;
	my @id=sort keys %{$hash_fa};
	for(my $i=0; $i<@id; $i++){
		my $line='';
		for(my $k=0;$k<length($hash_fa->{$id[$i]});$k++){
			my $l_tmp=substr($hash_fa->{$id[$i]},$k,1);
			if($l_tmp eq 'N'){if(special_Ns($id[$i],$k,1,$hash_fa)==1){$l_tmp=~tr/N/n/;}}
                	else{
                        	if(calculate_identity($id[$i],$k,1,$hash_fa)<$min_identity){$l_tmp=~tr/[ATCG]/[atcg]/;}
                	}
			$line.=$l_tmp;
		}
		while($line=~/([atcg]{$len_specific,})/g){my $L=length($1);my $a="n"x$L;$line=~s/$1/$a/;}
		$line=~tr/[atcg]/[ATCG]/;
		$hash_out1->{$id[$i]}=$line;
	}
}
######### get the low identity block. The window size is $windows. If meet '---' or 'n' or 'N', the window will be forward terminated.
sub low_identity{
	my($hash_out1,$hash_out2)=@_;
	my @id=sort keys %{$hash_out1};
	for(my $i=0;$i<@id;$i++){
		my $line=$hash_out1->{$id[$i]};
		my $line_tmp='';
		my @pos;$pos[0]=$pos[1]=$pos[2]=$pos[3]=0;##$pos[0] the start pos; $pos[1] the length; $pos[2] the flag 0:the prior is end by un'n', 1: is n.
		
		for(my $k=0;$k<length($line);$k++){
			my $base=substr($line,$k,1);
			$line_tmp.=$base;
			$pos[1]++;
			$pos[0]=$k if($pos[0]==0 && $pos[3]>0);
#			print "$pos[0]\n";
			if($base=~/n|N/){
				if($pos[2]==1 && $pos[1]<4){
					$line_tmp=~tr/[ATCG]/[atcg]/;
				}
				if(calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1)<$min_identity){
					$line_tmp=~tr/[ATCG]/[atcg]/;
				}
				$hash_out2->{$id[$i]}.=$line_tmp;
				$pos[0]=$pos[1]=0;$line_tmp='';$pos[2]=1;$pos[3]++;
				next;
			}
			
			if($line_tmp=~/-{6,}/ or $k+1==length($line)){
				if(calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1)<$min_identity){$line_tmp=~tr/[ATCG]/[atcg]/;}
				$hash_out2->{$id[$i]}.=$line_tmp;
				$pos[0]=$pos[1]=$pos[2]=0;$line_tmp='';$pos[3]++;
				next;
			}
			if(($pos[1]-$windows)==int(0.5*$windows+0.5)){
				$k=$k-int(0.5*$windows+0.5);
				$pos[1]=$pos[1]-int(0.5*$windows+0.5); $pos[3]++;
				my $flag_pos=$pos[1];
				my $li;
				my $identity_0=calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1);
				my $identity_tmp=$identity_0;
				if($identity_0<$min_identity){

					#while(1){
					#	$pos[1]--;$k--;
					#	my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1);
					#	if($identity_1>$identity_0){$identity_0=$identity_1;}
					#	else{last;}
					#}
					#$pos[1]++;$k++;
					my $pos_tmp=$pos[1];
					if($pos_tmp==$flag_pos){
						$identity_0=$identity_tmp;
						while(1){
							$pos[1]++;$k++;
							my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1);
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
						my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1);
						if($identity_1<$identity_0){$identity_0=$identity_1;}
						else{last;}
					}
					$pos[1]++;$k++;
					my $pos_tmp=$pos[1];
					if($pos_tmp==$flag_pos){
						$identity_0=$identity_tmp;
						while(1){
							$pos[1]++;$k++;
							my $identity_1=calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1);
							if($identity_1>$identity_0){$identity_0=$identity_1;}
							else{last;}
						}
						$pos[1]--;$k--;
					}
					$li=substr($line,$pos[0],$pos[1]);
				}
				$hash_out2->{$id[$i]}.=$li;
				$pos[0]=$pos[1]=$pos[2]=0;$line_tmp='';
				next;
			}
			
		}
		if($line_tmp){
			$pos[3]++;
			if(calculate_identity($id[$i],$pos[0],$pos[1],$hash_out1)<$min_identity){$line_tmp=~tr/[ATCG]/[atcg]/;}
			$hash_out2->{$id[$i]}.=$line_tmp;
		}
	}
}

sub special_Ns{
        my ($id,$start,$length,$hash_fa)=@_;
        my $f=1;
	my @id=sort keys %{$hash_fa};
        my $l1=substr($hash_fa->{$id},$start,$length);
        if($l1 ne 'N'){$f=0;}
        else{
                foreach(@id){
                        if($_ eq $id){next;}
                        my $l2=substr($hash_fa->{$_},$start,$length);
                        if($l2 ne '-'){$f=0;}
                }
        }
        return $f;
}

###compared with mini_identity
sub calculate_identity{
	my ($id,$start,$length,$hash_fa)=@_;
	my $f=0;
	my @id=sort keys %{$hash_fa};
	if($length==0){return $f;}
	my $l1=substr($hash_fa->{$id},$start,$length);
	foreach(@id){
		if($_ eq $id){next;}
		my $l2=substr($hash_fa->{$_},$start,$length);
		my $iden=get_identity($l1,$l2);
		if($iden>$f){$f=$iden;}
	}
	return $f;
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
	This program is used to mask the specific and low identity sequences of multiple sequence alignment in MAF format by 'n' and lowercase respectively.
	
	perl $0 --input <input.maf> --window <int number> --minidentity <float number> --len <int number>;
	or perl $0 -i <input.maf> -w <int number> -m <float number> -l <int number>;

	--input(i): multiple sequence alignment in MAF format;
	--window(w): int number,default 15;
	--minidentity(m):  minimum identity, default 0.3;
	--len(l): the minimum length of specific, default 6;
	--help(h): print the help information on screen.
	
	The output file: *.fliter:alignment file after masking;
	
USAGE
}
