'''
Created on Jun 30, 2009
First attempt at a data structure for high-resolution genome-wide data
@author: lgoff
'''
import genomelib
import gzip,time,sys
import copy
import numpy as np
from tables import *
import rpy
import Chip

class ContinuousData(object):
    '''
    Data storage object that is specific to a single chromosome
    '''
    
    def __init__(self,name,chr,binSize = 1,data = {}):
        '''
        Constructor: Creates instance specifically tailored to a given chromosome
        '''
        self.name = name
        self.chr = chr
        self.binSize = int(binSize)
        self.fname = "%s_%s_%d.bin" % (self.name,self.chr,self.binSize)
        if len(data.keys())>0:
            self.data = data
        else:
            self.data = {
                         '+':np.zeros(genomelib.chr_lengths[chr]/binSize,'d'),
                         '-':np.zeros(genomelib.chr_lengths[chr]/binSize ,'d')
                         }
    
    def __len__(self):
        """Equivalent to length of the genome"""
        return np.alen(self.data['+'])
    
    def __repr__(self):
        return self.name
    
    def __str__(self):
        return self.name
    
    def getMin(self,strand):
        return np.amin(self.data[strand])
    
    def getMax(self,strand):
        return np.amax(self.data[strand])
    
    def whichMax(self,strand):
        return np.argmax(self.data[strand])
    
    def whichMin(self,strand):
        return np.argmin(self.data[strand])
    
    def getDataRange(self,strand,start,end):
        return self.data[strand][(start/self.binSize)-1:(end/self.binSize)-1]
    
    def addInterval(self,interval):
        if self.chr != interval.chr:
            return "Wrong data file"
        else:
            self.data[interval.strand][(interval.start/self.binSize)-1:(interval.end/self.binSize)-1]=self.data[interval.strand][(interval.start/self.binSize)-1:(interval.end/self.binSize)-1]+interval.count
    
    def write(self,fname=None):
        if fname == None:
            fname = self.fname
        fd = gzip.open(fname,'wb')
        for s in self.data.keys():
            fd.write(self.data[s])
        fd.close()
    
    def read(self,fname):
        pass
    
    def innerHeight(self,strand,start,end):
        region = self.getDataRange(strand,start,end)
        return np.amax(region)
    
    def outerHeight(self,strand,start,end):
        region = self.getDataRange(strand,start,end)
        return sum(region)

class SimpleChIPData(object):
    
    def __init__(self,files):
        self.data = {}
        self.samples = []
        for fname in files:
            sampleName = fname.rstrip(".gff")
            self.samples.append(sampleName)
            sys.stderr.write("Parsing file '%s'...\n" % fname)
            self.data[sampleName] = Chip.parseNimblegen(fname)
        
    def doIt(self,permuted,windows=[5,6,7,8,9,10,11,12],threshold=0.05):
        self.normalize()
        self.joinProbes()
        for winSize in windows:
            self.scan(permuted,winSize,threshold)
    
    def makeMatrix(self):
        self.dataMatrix = np.empty((len(self.data[self.data.keys()[0]]),len(self.samples)),'f')
        for i in range(0,len(self.data.keys())):
            self.dataMatrix[:,i]=[x.score for x in self.data[self.data.keys()[i]]]
        sys.stderr.write("Created dataMatrix!\n")
    
    def quantileNormalize(self):
        if 'dataMatrix' not in self.__dict__: self.makeMatrix()
        rpy.r.library("limma")
        sys.stderr.write("Performing Quantile Normalization...\n")
        self.normMatrix = rpy.r.normalizeQuantiles(self.dataMatrix)
        
    def normalize(self):
        if 'normMatrix' not in self.__dict__: self.quantileNormalize()
        sys.stderr.write("Replacing values in data with normalized values...\n")
        for i in range(0,len(self.data.keys())):
            for j in range(0,np.shape(self.normMatrix)[0]):
                self.data[self.data.keys()[i]][j].score = self.normMatrix[j,i]
                
    def joinProbes(self):
        sys.stderr.write("Joining Probes into intervals...\n")
        self.intervals = {}
        for sample in self.samples:
            sys.stderr.write("\t%s\n" % sample)
            self.intervals[sample] = Chip.joinNimblegenIntervals(self.data[sample])
            
    def scan(self,permuted,windowSize,threshold=0.05):
        sys.stderr.write("Scanning with window of size %d..\n" % windowSize)
        for sample in self.samples:
            sys.stderr.write("\t%s\n" % sample)
            for i in self.intervals[sample]:
                i.scan(permuted,windowSize,threshold)
    
    
   