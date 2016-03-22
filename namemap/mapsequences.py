#! /usr/bin/env python
 
import os
import sys
import re

pattern = re.compile("(?<=[,(])([^,():]+)(?=[,():])")
 
if ("--help" in sys.argv) or ("-?" in sys.argv) or len(sys.argv) < 4:
    sys.stderr.write("usage: %s [<tree-file-path>] [<map-file-path>] [<out-file-path>] [-rev]\n	-rev: use for mapping from right hand side names in the mapping file (i.e. short names) to left hand side names (i.e. long names)\n"%sys.argv[0])
    sys.exit(1)
 
src_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
if not os.path.exists(src_fpath):
    sys.stderr.write('Not found: "%s"' % src_fpath)
src = open(src_fpath,"r")        
 
dest_fpath = os.path.expanduser(os.path.expandvars(sys.argv[3]))
dest = open(dest_fpath, "w")

map_fpath = os.path.expanduser(os.path.expandvars(sys.argv[2]))
mapfile = open(map_fpath, "r")

reverse = False
if len(sys.argv) == 5 and sys.argv[4] == "-rev":
    reverse = True

print "writing to file %s" %os.path.abspath(dest_fpath)

mapping = {}
for line in mapfile:
    m = line.replace("\n","").replace("\r","").split("\t")
    if not reverse:
        mapping[m[1].split(" ")[0]] = m[0]
    else: 
        mapping[m[0].split(" ")[0]] = m[1]
        
for t in src:    
    t = pattern.sub(lambda m: '"%s"' % mapping[m.group(1)],t)
    dest.write(t)

dest.close()
