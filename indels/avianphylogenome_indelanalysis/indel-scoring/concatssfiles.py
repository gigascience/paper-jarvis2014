#!/usr/bin/env python3
# Description: This script concatenates xread files, after trimming all but indel characters,
# into a single file that can be read with WinClada.
# Author: Gregory Penn gpenn@nmsu.edu
# Version: 0.9
# Date: 24 August 2012
# To do: Use temporary files for taxon files.

import argparse, re, glob, tempfile, sys, signal, os
# this handles the sigpipe signal sent by the shell if a pipe is interupted before writing is finished (as by head or tail)
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

parser = argparse.ArgumentParser(description="Merge ss files. Can size-select indels if ss files created by 2xread.")

parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.85")
parser.add_argument("src_path", metavar="path", type=str, help="Path to files to be merged; enclose in quotes, accepts * as wildcard for directories or filenames")
# parser.add_argument("-r", "--removenucleotides", action="store_true", help="Removes nucleotide characters, leaving only indels")
parser.add_argument("-l", "--log", action="store_true", help="Logs summary statistics")
parser.add_argument("-min", "--indellength_minimum", type=int, help="Minimum indel length. This option creates a matrix of indel characters only; it removes nucleotide characters.")
parser.add_argument("-max", "--indellength_maximum", type=int, help="Maximum indel length. This option creates a matrix of indel characters only; it removes nucleotide characters.")

args = parser.parse_args()

files = glob.glob(args.src_path)

if not files:
    print('File does not exist: ' + args.src_path, file=sys.stderr)

indelline_regexp = re.compile(r"^\{(\d+)\ssequence_indel_(\d+)-(\d+)")
# sequenceline = re.compile(r"^[A-Z]{5}\s+")
# sequenceline = re.compile(r"^([A-Z]{5})\s+([^\s]+)$")
sequenceline = re.compile(r"^([A-Z]{5})\s+([^\s]+)\s*$")
xreadline = re.compile(r"^xread\s*.*$")
matrixsizeline = re.compile(r"^([0-9]+)\s+([0-9]+)\s*$")

terminal_tmpfiles = {} # creates an empty dictionary

concatseq_length = 0
errorcount = 0
filecount = 0
noindels_count = 0
terminals_lengths = set([])

class Nucleotide:
    def __init__(self, file, terminal, position):
        pass

class Indel:
    def __init__(self, file, position, first, last, reverseindex):
        self.intron = file
        self.position = position
        self.first = first
        self.last = last
        self.length = last - first + 1
        self.reverseindex = reverseindex
        
class Terminal(str):
    def __init__(self, terminal):
        self.terminal = terminal
        self._introns = []
        

if args.indellength_minimum != None or args.indellength_maximum != None:
    if args.indellength_minimum != None and args.indellength_maximum == None:
        args.indellength_maximum = 10**9
    elif args.indellength_minimum == None and args.indellength_maximum != None:
        args.indellength_minimum = 0
    
    for file in files:
        foundmatrixsize = False
        indels = []
        
        with open(file, 'r') as thisfile:

            seekingmatrixsize = False

            for line in thisfile:
                if not foundmatrixsize and xreadline.match(line):   # Find matrix size
                    seekingmatrixsize = True
                if seekingmatrixsize and matrixsizeline.match(line):
                    sizeline = matrixsizeline.match(line)
                    seq_length = int(sizeline.group(1))
                    taxa_length = int(sizeline.group(2))
                    seekingmatrixsize = False
                    foundmatrixsize = True

                if indelline_regexp.match(line):   # Find indels' metadata
                    indelline = indelline_regexp.match(line)
                    position = int(indelline.group(1))
                    first = int(indelline.group(2)) - 1 # -1 required to correct stupid indexing.
                    last = int(indelline.group(3)) - 1  # ditto
                    length = last - first + 1
                    reverseindex = -1 * (seq_length - position)
                    indels.append(Indel(file, position, first, last, reverseindex))
        
        if not foundmatrixsize:
            err_msg = thisfile.name + ': No matrix found.'
            print(err_msg, file = sys.stderr)
        if foundmatrixsize == True:
            indellengthfilter = []
            if len(indels) > 0:
                for indel in indels:
                    if indel.length >= args.indellength_minimum and indel.length <= args.indellength_maximum:
                        indellengthfilter.append(indel.reverseindex)
                    thisseq_length = len(indellengthfilter)
                if thisseq_length == 0:
                    if args.log:
                        log_msg = thisfile.name + '\nNo indels meet criteria: minimum length = ' + str(args.indellength_minimum) + ', maximum length = ' + str(args.indellength_maximum) + '.\n'
                        with open('concat.log', 'a') as log:
                            log.write(log_msg)
                    
            else:
                thisseq_length = 0
                if args.log:
                    log_msg = thisfile.name + '\nIndel coding not found in file.'
                    with open('concat.log', 'a') as log:
                        log.write(log_msg)

            if thisseq_length > 0:
                with open(file, 'r') as thisfile:
        
                    # Note: reverse indexing gets around the problem of ambiguous characters,
                    # which have a length greater than 1 and through off normal forward indexing.
                    # This works because there are no ambiguous indel characters.
                    file_terminals = set([])
                    for line in thisfile:
                        if sequenceline.match(line):
                            terminal = sequenceline.match(line).group(1)
                            sequence = sequenceline.match(line).group(2)
                            file_terminals.add(terminal)
                            filteredindels = str()
                            for i in indellengthfilter:
                                filteredindels += sequence[i]
                            if terminal in terminal_tmpfiles: # check whether a tmpfile has been created for this terminal.
                                output = filteredindels
                            else:                    
                                tf = tempfile.NamedTemporaryFile(prefix=terminal) # create a tmpfile for this terminal.
                                terminal_tmpfiles[terminal] = tf.name # adds a key:value pair (terminal:tf.name) to the dictionary
                                tf.close()
                                prependdashes = "-" * concatseq_length # dashes get the new tmpfile up to length with others.
                                output = prependdashes + filteredindels

                            with open(terminal_tmpfiles[terminal], 'a') as tf:
                                tf.write(output)
                    for absentterminal in terminal_tmpfiles.keys() - file_terminals:
                        with open(terminal_tmpfiles[absentterminal], 'a') as tf:
                            output = '-' * thisseq_length
                            tf.write(output)
                    concatseq_length += thisseq_length
                    # write file and number of characters to log
                    # sys.stderr.write("hello\n")
                    if args.log:
                        with open('concat.log', 'a') as logfile:
                            logthisfile = str(file) + "\t" + str(concatseq_length - thisseq_length) + "\t" + str(concatseq_length - 1) + "\n"
                            logfile.write(logthisfile)
        else:
            sys.stderr.write("No matrix found in file: " + thisfile.name + '\n')
# if not removing nucleotide characters
else:
    for file in files:
        indels = []
        seq_length = None
        foundmatrixsize = False
        with open(file, 'r') as thisfile:

            seekingmatrixsize = False
            file_terminals = set([])

            for line in thisfile:
                if not foundmatrixsize and xreadline.match(line):   # Find matrix size
                    seekingmatrixsize = True
                if seekingmatrixsize and matrixsizeline.match(line):
                    sizeline = matrixsizeline.match(line)
                    seq_length = int(sizeline.group(1))
                    taxa_length = int(sizeline.group(2))
                    seekingmatrixsize = False
                    foundmatrixsize = True

                # Note: reverse indexing gets around the problem of ambiguous characters,
                # which have a length greater than 1 and through off normal forward indexing.
                # This works because there are no ambiguous indel characters.
                if sequenceline.match(line):
                    terminal = sequenceline.match(line).group(1)
                    sequence = sequenceline.match(line).group(2)
                    file_terminals.add(terminal)
                    if terminal in terminal_tmpfiles:
                        output = sequence
                    else:                    
                        tf = tempfile.NamedTemporaryFile(prefix=terminal)
                        terminal_tmpfiles[terminal] = tf.name
                        tf.close()
                        prependdashes = "-" * concatseq_length
                        output = prependdashes + sequence

                    with open(terminal_tmpfiles[terminal], 'a') as tf: # for this file
                        tf.write(output)
                        
        if foundmatrixsize:
            for absentterminal in terminal_tmpfiles.keys() - file_terminals:
                with open(terminal_tmpfiles[absentterminal], 'a') as tf:
                    output = '-' * seq_length
                    tf.write(output)
                    
            concatseq_length += seq_length
            
            # write index to log
            if args.log:
                with open('concat.log', 'a') as logfile:
                    logthisfile = str(file) + "\t" + str(concatseq_length - seq_length) + "\t" + str(concatseq_length - 1) + "\n"
                    logfile.write(logthisfile)
            
        else: sys.stderr.write("No matrix found in file: " + thisfile.name + '\n')

##########
if concatseq_length == 0:
    sys.exit("No output sent to stdout because concatenated matrix has length of 0.")
    
outlengths = set([])
outfile_terminals = set([])

line1 = 'xread\n'
line2 = str(concatseq_length) + ' ' + str(len(terminal_tmpfiles)) + '\n'
header = line1 + line2
endfile = ';\n'

sys.stdout.write(header)
    
for terminal in terminal_tmpfiles.keys():
    # print(terminal_tmpfiles[terminal])
    with open(terminal_tmpfiles[terminal], 'r') as src:
        srcseq = src.read()
        outlength = len(srcseq)
        outlengths.add(outlength)
        lineout = terminal + " " + srcseq + '\n'
        outfile_terminals.add(terminal)
        sys.stdout.write(lineout)
sys.stdout.write(endfile)
############
# if args.log:
#     log = str(filecount) + ' files\n'
#     log += str(len(taxon_tmpfiles)) + ' terminals\n'
#     if args.removenucleotides:
#         log += str(concatseq_length) + ' indel characters in matrix\n'
#         log += str(noindels_count) + ' files had no indels'
#     else:
#         log += str(concatseq_length) + ' characters in matrix'
#     with open('concatXread.log', 'w') as logfile:
#         logfile.write(log)
