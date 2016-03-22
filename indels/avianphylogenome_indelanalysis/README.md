# Avian Phylogenome Project - Indel Analysis #

This repository contains scripts and information pertaining to the indel analysis for the Avian Phylogenome Project. The analysis pipeline has very little source code and is dependent on various bioinformatics and systematics tools.

The analysis pipeline contains three distinct phases. For each phase, scripts and instructions are located within their respective directories:

1. 'indel-scoring' - Scoring indels/gaps (via Simple Indel Coding) in nucleotide alignments and generating indels alignment matrices filtered by indel size.
2. 'tree-estimation' - Estimating trees using indel alignments with parsimony and likelihood methods
3. 'node-RI' - Evaluating indel characters based on Retention Indices using trees generated with indel and nucleotide alignments.

### Dependencies ###

The scripts were developed and used within a Linux environment. These should work on Windows and Mac, as well, but I haven't tested them.

Before running the scripts, make sure the following packages are installed and available in your path:

* [2xread](http://www.nybg.org/files/scientists/2xread.html)
* [Python](https://www.python.org/) (Version 3)
* [TNT](http://www.lillo.org.ar/phylogeny/tnt/) (May 2013 or later)
* [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html)

### Comments, Bugs and Suggestions ###

I created these scripts as I learned scripting and programming, and to serve the immediate needs of the project. Hence, they are not as portable as they could be. Please contact me at nitishnarula19@gmail.com.  For the latest version of these scripts, see: https://nitishnarula@bitbucket.org/nitishnarula/avianphylogenome_indelanalysis.git
