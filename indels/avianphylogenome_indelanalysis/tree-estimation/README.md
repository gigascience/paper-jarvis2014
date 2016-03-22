# Tree Estimation with Indel Matrices #

### Parsimony Analysis with TNT ###

[TNT](http://www.lillo.org.ar/phylogeny/tnt/) was be used for maximum parsimony tree estimation. It can handle large alignment matrices like those generated in this project.

The fulltnt.run file is a TNT script file used in this project. The script runs 1000 multiple random addition searches, saves the consensus of the best trees, runs 10000 bootstraps
and saves the bootstrap support values on the consensus tree. You can run the script in the TNT shell, as follows:

	$ tnt
	tnt*>mxram 2000;
	tnt*>fulltnt indelmatrix.ss;

Notes:

- You **must** provide suitable ram to TNT to open your alignment.

- The alignment file (or matrix file) **must** contain 'proc/;' (without the single quotes) at the end of the file. For our pipeline, indel matrices generated in the previous step were appended with this token as follows:

		$ echo "proc/;" >> indelmatrix.ss

- stats.run script, provided by TNT authors, **must** be present in the same directory as the fulltnt.run script.

- Change the output file names (mynelsen.tre, finalconsensus.tre and logfile.txt) in the TNT script according to your needs.

For more on TNT scripting, see the TNT scripts [page](http://www.lillo.org.ar/phylogeny/tnt/scripts/).

TNT trees are not well suited for most tree visualization packages. Use the tnt2newick script available on [TNT Wiki](http://tnt.insectmuseum.org/index.php/FAQ#How_can_I_export_trees_in_TNT_format_then_convert_them_to_Newick_format_in_an_automated_way.2C_with_a_Python_script.3F)

### Likelihood Analysis with RAxML ###

[RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html) was used for maximum likelihood tree estimation. Like TNT, it is capable of handling large datasets.

The following three steps were carried out for the ML trees generated in the indel analysis. The first two can run separately at the same time. The last step is run after the first two are complete:

* Likelihood searches

		$ raxml -T [number of cores] -m BINGAMMA -p $RANDOM -N [number of searches, typically 10] -s [input phylip matrix] -w [output directory] -n [output name]

* Bootstrapping

		$ raxml -T [number of cores] -m BINGAMMA -b $RANDOM -p $RANDOM -N autoMRE -s [input phylip matrix] -w [output directory] -n [output name]

* Best Tree with bootstrap support

		$ raxml -T [number of cores] -m BINGAMMA -f b -t [RAxML bestTree from step 1] -z [RAxML bootstrap tree from step2] -w [output directory] -n [output name]
