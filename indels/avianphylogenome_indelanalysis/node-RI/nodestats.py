#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  nodestats.py
#  
#  Copyright 2013 Nitish Narula <nnarula@nmsu.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import os, sys, re

if ("--help" in sys.argv) or ("-?" in sys.argv) or len(sys.argv) < 2:
    sys.stderr.write("usage: %s [<charRI output file-path>]\n"%sys.argv[0])
    sys.exit(1)    
 
src_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
if not os.path.exists(src_fpath):
    sys.stderr.write('Not found: "%s"' % src_fpath)
src = open(src_fpath, "r")


def nodestats(RIlist):
	if len(RIlist) == 0:
		return "0, 0.0, 0.0"
	else:
		counter = 0
		RIsum = 0.000
		for i in range(0, len(RIlist)):
			if float(RIlist[i]) == 1.000:
				counter += 1
			RIsum += float(RIlist[i])
		result = str(len(RIlist)) + ", " + str(RIsum/len(RIlist)) + ", " + str(float(counter)/float(len(RIlist)))
		return result


searchtree = re.compile(r"Tree\s\d+\s:")
termtaxon= re.compile(r"\s+[A-Z]{5}\s:")
nodenum = re.compile(r"\s+Node\s\d+\s:")
charline = re.compile(r"\s+Char.\s\d+:")
charRI = re.compile(r"\d{1}\.\d{3}")
RIlist = []
readRI = False

#char_list = char_lines.findall(src)
for line in src:
	if searchtree.match(line):
		print(line, end=' ') #searchtree.match(line).group()
	elif termtaxon.match(line) or nodenum.match(line):
		if not readRI:
			readRI = True
		else:
			print(thisnode, nodestats(RIlist))
			RIlist[:] = []
		if termtaxon.match(line):
			thisnode = termtaxon.match(line).group()
		else:
			thisnode = nodenum.match(line).group()
	elif readRI:
		if charline.match(line):
			RIlist.append(charRI.search(line).group())

print(thisnode, nodestats(RIlist))
src.close()
