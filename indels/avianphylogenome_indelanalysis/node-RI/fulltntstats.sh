#!/bin/bash
# $Id$

usage ()
{
    (
	echo "${0} -- Run stats and aporis on matrices given a tree using TNT."
	echo "usage: ${0} [-h] -i <xread file> -t <tree file> -o <output directory> -n <name without spaces>"
    	echo "TNT executable must be in your path. TNT scripts - aporis.run and stats.run - must be in the current directory."
    	echo "Change 'mxram' and 'txtsize' arguments in this script to give more RAM and output size to TNT."
    ) > /dev/stderr
    exit 1
}

cleanup()
{
    rm -f $RUNFILE
}

trap cleanup HUP INT QUIT TERM 

: ${TMPDIR:=/tmp}

RUNFILE=${TMPDIR}/$$.run

# command line arguments
INFILE=
OUTDIR=
TREEFILE=
PROJNAME=

while getopts "hi:t:o:n:" opt; do
  case $opt in
    h) usage;;
    i) INFILE=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) TREEFILE=$OPTARG;;
    n) PROJNAME=$OPTARG;;
    \?) usage;;
  esac
done 

if [ "X$INFILE" = "X" || "X$OUTDIR" = "X" || "X$TREEFILE" = "X" ]; then 
    echo "Missing argument(s)" > /dev/stderr
    usage
fi

$CAT << EOF > $RUNFILE
mxram 4000;
silent = buffer;
proc $INFILE;
xinact;
xread -;
proc $TREEFILE;
txtsize 20480;
apo -;
log $OUTDIR/$PROJNAME.apo.txt;
svtxt;
log/;
clbuffer;
aporis 0 $OUTDIR/$PROJNAME.aporis.txt;
log $OUTDIR/$PROJNAME.stats.txt;
stats;
svtxt;
log/;
clbuffer;
taxname=;
tsave *$OUTDIR/$PROJNAME.blen.tre;
ttags=;
blen *0;
save *;
tsave/;
z;
EOF

tnt p $RUNFILE

cleanup
