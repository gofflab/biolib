'''
Created on Mar 17, 2011

@author: lgoff
'''
# from misc import rstrips  # rasmus library removed - not Python 3.12 compatible

def bed2GTF(fname,outfile=None):
    """This does not work yet"""
    handle = open(fname,'r')
    if outfile == None:
        outfile = fname.rstrip('.bed')+'.gtf'
    outHandle = open(outfile,'w')
    for line in handle:
        line = line.rstrip()
        if line.startswith("#"):
            print(line, file=outHandle)
            continue
        if line.startswith("track") or line.startswith("browser"):
            print(line, file=outHandle)
            continue
        vals = line.split("\t")
    pass
