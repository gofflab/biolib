'''
File format conversion utilities for genomic annotation files.

Contains functions for converting between common bioinformatics file formats
such as BED and GTF.

Created on Mar 17, 2011

@author: lgoff
'''
# from misc import rstrips  # rasmus library removed - not Python 3.12 compatible

def bed2GTF(fname, outfile=None):
    """Convert a BED file to GTF format (not yet fully implemented).

    Opens the input BED file, writes comment lines and track/browser header
    lines through unchanged, and parses remaining tab-delimited lines. The
    actual record conversion logic is not yet implemented.

    Note: This function is incomplete and does not currently produce GTF
    output records.

    Args:
        fname: Path to the input BED file.
        outfile: Path for the output GTF file. Defaults to fname with the
            trailing '.bed' stripped and '.gtf' appended.
    """
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
