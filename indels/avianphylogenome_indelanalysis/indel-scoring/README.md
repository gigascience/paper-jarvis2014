# Indel Scoring and Matrices #

### Scoring with 2xread ###
Indels, or gaps, are scored with the Perl script [2xread](http://www.nybg.org/files/scientists/2xread.html). Make sure you have it available in your path.

The script can only score FASTA or dread alignment files. It generates an xread-formatted output (or Hennig86 fomrat) to standard out. This format can be read by packages like TNT, NONA and WinClada.

Here is a sample command to be run on the terminal:

	$ perl 2xread my_fasta_aln.fasta > output_file.ss

Note: .ss extension is the typical method for declaring a file as xread formatted.

2xread can take a second optional argument which the author it calls stem. This is used to define what kind of alignment it is - DNA or Amino Acid, etc. 
The defualt is the string "sequence". This option does not change how the script works; it is just for scripting purposes.

If there are many fasta files to be scored, a shell 'for' loop can come in handy. However, you will have to deal with output file names. Here's one way to do it:

	$ for fasta in fasta_dir/*.fasta
	do
		token = $(echo `basename $fasta | sed 's/.fasta//')
		perl 2xread $fasta > $token"_indels.ss"
	done 

### Concatenation and Size Filtration with concatssfiles.py ###

The Python3 script, concatssfiles.py, included here, is a custom script created by Gregory Penn to help concatenate many xread files.
Moreover, the script can filter out indel characters based on size.

Note that this script was written specifically for the Avian Phylogenome Project. 
For instance, in the project taxon names were standardized to 5 uppercase letters so as to easily access sequence identifiers programmatically.
This script utilizes this to identify lines containing sequences.

Call the script with '--help' option to see how it works.

### Converting to phylip format ###

For many phylogenetic applications, the xread format is not suitable. RAxML, for instance, requires the commonly used phylip format. 

To convert xread to phylip, the following 'awk' script can be:

	$ awk 'NR==2 {printf " %s %s\n", $2, $1} NR > 2 {print $0}' indel_matrix.ss > indel_matrix.phy