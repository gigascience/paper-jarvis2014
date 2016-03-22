#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  charRI.py
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

if ("--help" in sys.argv) or ("-?" in sys.argv) or len(sys.argv) < 3:
    sys.stderr.write("usage: %s [<apo output filepath>] [<aporis output filepath>]\n"%sys.argv[0])
    sys.exit(1)    
 
src_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
if not os.path.exists(src_fpath):
    sys.stderr.write('Not found: "%s"' % src_fpath)
src = open(src_fpath, "r")

ri_fpath = os.path.expanduser(os.path.expandvars(sys.argv[2]))


char_lines = re.compile(r"^\s+Char.*")
just_char = re.compile(r"Char.\s\d+:")
ri_str = re.compile(r"RI\s\d\.\d+")

#char_list = char_lines.findall(src)
for line in src:
	if  char_lines.match(line) is None:
		print(line, end=' ')
	else:
		the_line = char_lines.match(line).group()
		#if isinstance(the_line, basestring):
			#print "yes"
		tosearch = just_char.search(the_line).group()
		with open(ri_fpath, "r") as RIfile:
			for other_line in RIfile:
				if tosearch in other_line:
					ri_print = ri_str.search(other_line).group()
					print("     ", tosearch, ri_print)



src.close()
