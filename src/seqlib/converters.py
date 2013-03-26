'''
Created on Mar 17, 2011

@author: lgoff
'''
from misc import rstrips

def bed2GTF(fname,outfile=None):
    """This does not work yet"""
    handle = open(fname,'r')
    if outfile == None:
        outfile = rstrips(fname,'.bed')+'.gtf'
    outHandle = open(outfile,'w')
    for line in handle:
        line = line.rstrip()
        if line.startswith("#"):
            print >>outHandle, line
            continue
        if line.startswith("track") or line.startswith("browser"):
            print >>outHandle, line
            continue
        vals = line.split("\t")
    pass