'''
Created on Jul 6, 2009
This module will attempt to deal with the nimblegen array data in a similar mechanism to that achieved by Guttman et al.

@author: lgoff
'''
import Alignment,copy,rpy,random
import numpy as np
from intervallib import *
from misc import pp
import sys,glob
import continuousData

class ChipInterval(Interval):
    """Extends basic Interval class with Tiling array methods and attributes"""
    
    def __init__(self, chr, start, end, strand="*", score=0.0, readcount = -1,name="",sequence = "",data={}):
        Interval.__init__(self, chr, start, end, strand=strand, score=score, readcount = readcount,name=name,sequence = sequence,data=data)
        self.parents = []
        self.children = []
        
    def addChild(self, child):
        """Adds child node to self.children"""
        #assert child not in self.children
        if child not in self.children:
            child.parents.append(self)
            self.children.append(child)
    
    def removeChild(self, child):
        """Removes child node from self.children (not sure how or if this works. Don't trust it yet)"""
        child.parents.remove(self)
        self.children.remove(child)
    
    def childScores(self):
        """Returns list of scores for each interval in self.children"""
        return [x.score for x in self.children]
    
    def childAvg(self):
        """Empty"""
        pass
    
    def childMedian(self):
        """Empty"""
        pass
    
    def makeValMap(self,value = 'readcount'):
        """Check these two to see which one is right..."""
        self.valMap = np.zeros(len(self))
        self.valMap = self.valMap-1
        myTmp = []
        for x in range(0,len(self)):
            myTmp.append([])
        for i in self.children:
            for j in range(i.start,i.end+1):
                myTmp[j-self.start].append(i.__dict__[value])
        for nt in range(0,len(myTmp)):
            if len(myTmp[nt])>0:
                self.valMap[nt]=sum(myTmp[nt])/len(myTmp[nt])

    
    """
    #This does not work at all....    
    def makeValMap(self):
        
        self.valMap = np.zeros(len(self))
        self.valMap = self.valMap-1
        for i in self.children:
            for j in range(i.start,i.end+1):
                if self.valMap[j-self.start]<0:
                    self.valMap[j-self.start]=i.score
                else:
                     self.valMap[j-self.start]=(self.valMap[j-self.start]+i.score)/2
    
    
    def makeValMap(self):
        '''Check these two to see which one is right...'''
        self.valMap = np.zeros(len(self))
        self.valMap = self.valMap-1
        myTmp = []
        for x in range(0,len(self)):
            myTmp.append([])
        for i in self.children:
            for j in range(i.start,i.end+1):
                myTmp[j-self.start].append(i.score)
        for nt in range(0,len(myTmp)):
            if len(myTmp[nt])>0:
                self.valMap[nt]=sum(myTmp[nt])/len(myTmp[nt])
        #pp(myTmp,1)        
    """
    
    def plotVals(self):
        """Creates a line plot (via rpy) across all bases within interval of the scores from self.valMap for the given base"""        
        if 'valMap' not in self.__dict__:
            self.makeValMap()
        rpy.r.x11()
        #rpy.r.plot(range(self.start,self.end+1),self.valMap,ylab="",type="l",lwd=2,main=str(self))
        rpy.r.plot((self.children[0].start,self.children[0].end),(self.children[0].score,self.children[0].score),type="l",lwd = 2,ylim=(min(c.score for c in self.children),max(c.score for c in self.children)))
        for x in self.children[1:]:
            rpy.r.lines((x.start,x.end),(x.score,x.score),lwd=2)
        
    def plot(self):
        """Convenience wrapper for self.plotVals"""
        self.plotVals()
    
#    def uniqifySig(self):
#        keys = {}
#        for e in self.significant:
#            keys[e] = 1
#        self.significant = keys.keys()
    
    def scan(self,permuted,windowSize,threshold):
        self.children.sort()
        if 'significant' not in self.__dict__: 
            self.significant = []
        for i in range(0,len(self.children)-windowSize):
            tester = np.mean([x.score for x in self.children[i:i+windowSize]])
            if len(permuted[windowSize][permuted[windowSize]>=tester])/float(len(permuted[windowSize]))<=threshold:
                for j in self.children[i:i+windowSize]:
                    if j not in self.significant:
                        k = copy.copy(j)
                        k.children = []
                        self.significant.extend(j)
        
        

#This should be deleted...
class ChipData(object):
    """Container for one array's worth of NimbleGen data"""
    def __init__(self, fname, sampleName):
        self.fname = fname
        self.sampleName = sampleName
        self.probeData = {}
        
        #Populate self.probeData
        ChipIter = parseNimblegen(fname)
        for ci in ChipIter:
            if not ci.chr in self.probeData.keys():
                self.probeData[ci.chr] = []
            self.probeData[ci.chr].append(ci)
    
    def sort(self):
        """Sorts all chromosomes seperately and in place"""
        for k in self.data.keys():
            self.data[k].sort()
    
    def shuffle(self,chr):
        """This doesn't work yet"""
        vals = [x.score for x in self.probeData[chr]]
        return random.shuffle(vals)
            
#End crap   
       
def nimblegenIter(fname):
    """Returns a generator of ChipInterval objects from a nimblegen .GFF output file"""
    handle = open(fname,'r')
    for line in handle:
        if line.startswith("#"): continue
        tokens = line.split("\t")
        pname = tokens[8].split(";")[1].split("=")[1]
        yield ChipInterval(tokens[0],tokens[3],tokens[4],score=tokens[5],name=pname)
    
def parseNimblegen(fname):
    iter = nimblegenIter(fname)
    rtrn = []
    for i in iter:
        rtrn.append(i)
    return rtrn

def joinNimblegenIntervals(intervals,start='start',end='end',offset=1000):
    """
    Returns a list of independent transcription units overlaping by offset
    """
    
    if not intervals: return intervals
    
    intervals.sort()
    
    non_overlapping = []  
    current = copy.copy(intervals[0])
    current.addChild(copy.copy(current))
    current.score = 0.0
    for x in intervals[1:]:
        next = copy.copy(x)
        if current.intersects(next,start = start,end = end, offset = offset):
            current.end = max(current.end,next.end)
            current.addChild(copy.copy(next))
        else:
            current.name=""
            non_overlapping.append(current)
            current = copy.copy(next)
            current.addChild(copy.copy(current))
            current.score = 0.0
    current.name = ""
    non_overlapping.append(current)
    return non_overlapping

def probeScores(probes):
    """Returns list of scores across all a list of probes"""
    return np.array([x.score for x in probes],dtype='f')

def getRandomDist(probes,nRandom,windowSize):
    """Returns a numpy array of length 'nRandom' corresponding to the max values of sliding windows of size 'windowSize'
    from shuffled probe data.
    """
    sys.stderr.write("Getting %d Max value distributions from windows of size %d:\n" % (nRandom,windowSize))
    #scores = probeScores(probes)
    scores = probes
    lenScores = len(scores)
    maxVals = np.empty(nRandom,dtype = 'f')
    sys.stderr.write("\tIterations:\n")
    for i in range(0,nRandom):
        if i%10 == 0: sys.stderr.write("\t%d\n" % i)
        random.shuffle(scores)
        myMax = 0.0
        for win in range(0,lenScores-windowSize):
            avg = np.mean(scores[win:win+windowSize])
            if avg > myMax: myMax = avg
        maxVals[i] = myMax
    return maxVals

def calcPVals(segScores,permuted,windowSize):
    """This does not work yet"""
    return len(permuted[permuted>=segScores])


def main():
    files = glob.glob("*.gff")
    data = continuousData.SimpleChIPData(files)
    data.normalize()
    data.joinProbes()
    permuted = {}
    windows = [5,7,9,11]
    sys.stderr.write("Generating random data permutations...\n")
    for windowSize in windows:
        sys.stderr.write("\t%d\n" % windowSize)
        permuted[windowSize] = getRandomDist(data.data[data.samples[0]],1000,windowSize)
    
if __name__=="__main__":
    main()
            