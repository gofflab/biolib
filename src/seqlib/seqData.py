#!/usr/bin/env python
'''
Created on Oct 27, 2009

@author: lgoff
'''
import pysam
import mySam
from rpy import *
import copy
import intervallib

class SamData:
    def __init__(self,name,file,description):
        self.name = name
        self.file = file
        self.description = description
        self.type = 'basic'
        self.open()
    
    def __str__(self):
        return self.name
    
    def open(self):
        """Returns a pysam handle to the .BAM file"""
        self.handle = pysam.Samfile(self.file,'rb')
    
    def close(self):
        self.handle.close()
    
    def samSort(self):
        pass
    
    def samIndex(self):
        pass

    def pileupQuery(self,chr,start='',end=''):
        pos = []
        n = []
        for pileupcolumn in self.handle.pileup(chr,start,end):
            pos.append(pileupcolumn.pos)
            n.append(pileupcolumn.n)
        return (pos,n)

class ChromData(SamData):
    def __init__(self,name,file,description,mark,cellLine):
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
    """Incorporates strandedness and possibly an extension factor to account for fragment size"""
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
    """Not terribly flexible at this point, but will plot 'tracks' from a given chrom,start,end 
    position from a list of opened .BAM files"""
    
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
    '''
    get the bitvalue of denary (base 10) number n at the equivalent binary
    position p (binary count starts at position 0 from the right)
    '''
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
    strand = strandFlag(int(samRead.flag))
    return intervallib.Interval(samRead.qname,int(samRead.pos)+1,int(samRead.pos)+samRead.rlen+1,strand)

def samReads2Intervals(samReads,start='start',end='end',score='readcount',sampleName=".",offset=0):
    """samReads is an iterator object over a set of sam reads using the pysam 'fetch' call"""
    pass

    