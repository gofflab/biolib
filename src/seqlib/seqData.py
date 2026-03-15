#!/usr/bin/env python
"""Data structures and utilities for working with BAM/SAM sequencing data.

Provides SamData and ChromData classes wrapping pysam for read access to BAM
files, a plotRegions function for strand-aware coverage visualisation via rpy,
and helper utilities for parsing SAM bitflags and converting reads to Interval
objects.

Note: This module depends on pysam and rpy, which must be installed separately.

Originally created on Oct 27, 2009.

Author: lgoff
"""

import intervallib
import pysam
from rpy import *


class SamData:
    """Wrapper around a pysam BAM file handle.

    Provides basic access to a sorted, indexed BAM file including pileup
    queries and a pysam Samfile handle.

    Attributes:
        name: Sample name string.
        file: Path to the BAM file.
        description: Human-readable description string.
        type: Data type label (default "basic").
        handle: Open pysam.Samfile handle.
    """
    def __init__(self,name,file,description):
        """Initialize and open a SamData object.

        Args:
            name: Sample name string.
            file: Path to the BAM file.
            description: Human-readable description of the sample.
        """
        self.name = name
        self.file = file
        self.description = description
        self.type = 'basic'
        self.open()

    def __str__(self):
        """Return the sample name string."""
        return self.name

    def open(self):
        """Open the BAM file and store the pysam handle in self.handle."""
        self.handle = pysam.Samfile(self.file,'rb')

    def close(self):
        """Close the pysam BAM file handle."""
        self.handle.close()

    def samSort(self):
        """Placeholder for BAM sorting (not yet implemented)."""
        pass

    def samIndex(self):
        """Placeholder for BAM indexing (not yet implemented)."""
        pass

    def pileupQuery(self,chr,start='',end=''):
        """Return per-position pileup depths for a genomic region.

        Args:
            chr: Chromosome name string.
            start: Start coordinate (default "" for beginning of chromosome).
            end: End coordinate (default "" for end of chromosome).

        Returns:
            A tuple (pos, n) where pos is a list of genomic positions and
            n is a list of corresponding pileup depths.
        """
        pos = []
        n = []
        for pileupcolumn in self.handle.pileup(chr,start,end):
            pos.append(pileupcolumn.pos)
            n.append(pileupcolumn.n)
        return (pos,n)

class ChromData(SamData):
    """SamData subclass for chromatin modification ChIP-seq BAM files.

    Extends SamData with mark and cell-line metadata.

    Attributes:
        mark: Histone mark or chromatin feature name (e.g. "H3K4me3").
        cellLine: Cell line identifier string.
        type: Data type label (always "chromatin").
    """
    def __init__(self,name,file,description,mark,cellLine):
        """Initialize a ChromData object.

        Args:
            name: Sample name string.
            file: Path to the BAM file.
            description: Human-readable description.
            mark: Histone mark or antibody target name.
            cellLine: Cell line identifier string.
        """
        SamData.__init__(self, name=name, file=file, description=description)
        self.mark = mark
        self.cellLine = cellLine
        self.type = 'chromatin'

##########
#Chromatin Data
##########

#UCSD Reference Epigenome (GSE16256, H1, various modifications)
ucsd_basedir = '/seq/compbio-hp/lgoff/raw_sequencing_data/chromatin/UCSC_Reference_Epigenome/GSE16256/'
GSE16256 = {'H3K4me1': ucsd_basedir+'GSM409307_UCSD.H3K4me1_sorted.bam',
            'H3K9me3': ucsd_basedir+'GSM409310_UCSD.H3K9me3_sorted.bam',
            'H3K27me3': ucsd_basedir+'GSM434776_UCSD.H3K27me3_sorted.bam',
            'H3K4me3': ucsd_basedir+'GSM409308_UCSD.H3K4me3_sorted.bam',
            'H3K36me3': ucsd_basedir+'GSM409312_UCSD.H3K36me3_sorted.bam',
            'H3K9ac': ucsd_basedir+'GSM434785_UCSD.H3K9ac_sorted.bam'
            }

#Broad Reference Epigenome (GSE17312)
broad_ref_basedir = '/seq/compbio-hp/lgoff/raw_sequencing_data/chromatin/Broad_Reference_Epigenome/GSE17312/'
GSE17312 = {'H3K27me3': broad_ref_basedir+'GSM433167_BI.H3K27me3_sorted.bam',
            'H3K9ac': broad_ref_basedir+'GSM433171_BI.H3K9ac_sorted.bam',
            'H3K36me3': broad_ref_basedir+'GSM433176_BI.H3K36me3_sorted.bam',
            'H3K4me3': broad_ref_basedir+'GSM433170_BI.H3K4me3_sorted.bam',
            'H3K9me3': broad_ref_basedir+'GSM433174_BI.H3K9me3_sorted.bam',
            'H3K4me1': broad_ref_basedir+'GSM433177_BI.H3K4me1_sorted.bam'
            }

def openBams(dataDict,cellLine):
    """Open a collection of BAM files described by a dictionary.

    Creates ChromData objects for each entry in dataDict, opens each BAM
    file handle, and returns the list.

    Args:
        dataDict: Dict mapping mark name to BAM file path.
        cellLine: Cell line identifier assigned to all ChromData objects.

    Returns:
        List of opened ChromData objects.
    """
    files = []
    for k,v in dataDict.items():
        sample = v.split("_")[0]
        files.append(ChromData(k,v,sample,k,cellLine))
    for f in files:
        f.open()
    return files

###############
#Utilities for dealing with sorted and indexed .BAM files
##############

"""

def plotRegions(bamHandle,chrom,start,end):
    '''Depricated'''
    pos = []
    n = []
    for column in bamHandle.pileup(chrom,start,end):
        pos.append(column.pos)
        n.append(column.n)
    r.plot(pos,n,type="h",col="purple",xlab=chrom+" position", ylab = "Aligned Reads", xlim=[start,end],main="Reads")

"""
def plotRegions(bamHandle,chrom,start,end):
    """Plot strand-aware read coverage for a genomic region using rpy.

    Counts per-position forward ("+") and reverse ("-") read coverage using
    pysam fetch, then draws a coverage plot via rpy with forward reads in blue
    above the axis and reverse reads in red below.

    Args:
        bamHandle: An open pysam Samfile or AlignmentFile handle.
        chrom: Chromosome name string.
        start: Start coordinate (integer).
        end: End coordinate (integer).
    """
    tmp = {}
    tmp["+"] = {}
    tmp["-"] = {}
    for read in bamHandle.fetch(chrom,start,end):
        if read.is_reverse == 0:
            for i in range(read.pos+1,read.pos+1+len(read.seq)):
                tmp["+"][i] = 1 + tmp["+"].get(i,0)
        elif read.is_reverse == 1:
            for i in range(read.pos+1,read.pos+1+len(read.seq)):
                tmp["-"][i] = 1 + tmp["-"].get(i,0)
    try: max_cov = max(tmp['+'].values()+tmp['-'].values())
    except ValueError: max_cov = 1

    r.plot(tmp['+'].keys(),tmp['+'].values(),type="h",col = "blue", ylim=[-max_cov,max_cov], xlab = chrom+" position", ylab = "Align Reads", xlim=[start,end], main = "Coverage "+chrom+":"+str(start)+"-"+str(end))
    r.lines(tmp['-'].keys(),map(lambda x: -x,tmp['-'].values()),type="h",col="red")
    r.abline(h=0,col="grey")


def plotChromProfile(bamFiles,chrom,start,end):
    """Plot stacked pileup-depth tracks for multiple BAM files via rpy.

    Opens a new rpy graphics device and plots one coverage track per BAM
    file in a vertically stacked layout. Not very flexible at this point.

    Args:
        bamFiles: List of opened SamData (or similar) objects with a
            .handle attribute supporting pileup() and a .name attribute.
        chrom: Chromosome name string.
        start: Start coordinate (integer).
        end: End coordinate (integer).
    """

    r.x11(width=6,height=10)
    r.par(mfrow=[len(bamFiles),1])
    for fname in bamFiles:
        pos = []
        n = []
        for column in fname.handle.pileup(chrom,start,end):
            pos.append(column.pos)
            n.append(column.n)
        r.plot(pos,n,type="h",xlab=chrom+" position",ylab="Aligned Reads",xlim=[start,end],ylim=[0,12],main=fname.name)

###############
#Functions for sam Reads
###############
def getBitValue(n, p):
    """Return the bit value of integer n at binary position p.

    Binary position 0 is the least significant bit (rightmost).

    Args:
        n: Denary (base-10) integer.
        p: Bit position to inspect (0-indexed from the right).

    Returns:
        0 or 1 depending on the bit at position p.
    """
    return (n >> p) & 1

def strandFlag(flag):
    """Returns strand of sequence from SAM record bitflag (field 4)"""
    flag = int(flag)
    if getBitValue(flag,4)==0:
        return "+"
    elif getBitValue(flag,4)==1:
        return "-"
    else:
        return "*"

def samRead2Interval(samRead):
    """Convert a single pysam AlignedRead to an intervallib.Interval.

    The strand is determined from the SAM bitflag. Coordinates are converted
    to 1-based by adding 1 to samRead.pos.

    Args:
        samRead: A pysam AlignedRead object.

    Returns:
        An intervallib.Interval with chr set to samRead.qname, 1-based
        start/end coordinates, and strand derived from the bitflag.
    """
    strand = strandFlag(int(samRead.flag))
    return intervallib.Interval(samRead.qname,int(samRead.pos)+1,int(samRead.pos)+samRead.rlen+1,strand)

def samReads2Intervals(samReads,start='start',end='end',score='readcount',sampleName=".",offset=0):
    """Convert a pysam fetch iterator of SAM reads to Interval objects.

    Note: This function is not yet implemented (passes without action).

    Args:
        samReads: Iterator object over SAM reads from a pysam 'fetch' call.
        start: Name of the start coordinate field (default "start").
        end: Name of the end coordinate field (default "end").
        score: Name of the score field (default "readcount").
        sampleName: Sample name string (default ".").
        offset: Integer offset applied to coordinates (default 0).
    """
    pass

