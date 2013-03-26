'''
Created on Jun 30, 2009

@author: lgoff
'''
from intervallib import *
import misc

class Alignment(object):
    """
    Basic Alignment class for short RNA reads
    Can be avoided directly in favor of aligner-specific implementations (ie. ShrimpRead and/or MAQRead)
    """
    def __init__(self,readname,chr,start,end,strand,score=0,readcount = -1,readsequence=''):
        self.readname = str(readname)
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.score = float(score)
        self.readsequence = readsequence
        self.readcount = readcount
    
    def __cmp__(self,b):
        return -cmp(self.score,b.score)    
    
    def __str__(self):
        return "%s:%s:%d:%d" % (self.readname,self.chr,self.start,self.end)
    
    def __repr__(self):
        return "%s:%s:%d:%d" % (self.readname,self.chr,self.start,self.end)
    
    def __len__(self):
        return self.end-self.start+1
    
    def isPlus(self):
        if self.strand=="+":
            return True
        else:
            return False
        
    def isMinus(self):
        if self.strand=="-":
            return True
        else:
            return False
    
    def toInterval(self):
        return Interval(self.chr,self.start,self.end,self.strand,self.score,self.readcount,name=self.readname)
    
    def toBed(self):
        return ("%s\t%d\t%d\t%s\t%d\t%s\n" % (self.chr,self.start,self.end,misc.seq2nuID(self.readsequence),self.readcount,self.strand))