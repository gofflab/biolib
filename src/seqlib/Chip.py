'''
Tools for working with NimbleGen ChIP-chip tiling array data.

Provides an interval class with tiling-array-specific methods, parsers for
NimbleGen GFF output files, interval-merging utilities, and statistical
helpers for identifying enriched regions via permutation-based p-value
estimation — following the approach of Guttman et al.

Created on Jul 6, 2009

@author: lgoff
'''
import copy
import glob
import random

# from misc import pp  # rasmus library removed - not Python 3.12 compatible
import sys

import numpy as np
import rpy2.robjects as robjects

from . import continuousData
from .intervallib import *


class ChipInterval(Interval):
    """Genomic interval extended with tiling-array probe-hierarchy support.

    Extends the basic Interval class with parent/child relationships so
    that individual NimbleGen probes (children) can be grouped under a
    merged enriched region (parent), and provides methods for computing
    coverage maps and plots from the probe scores.

    Attributes:
        parents: List of ChipInterval objects that contain this interval.
        children: List of ChipInterval objects contained within this
            interval (e.g. individual probes belonging to an enriched
            region).
    """

    def __init__(self, chr, start, end, strand="*", score=0.0, readcount = -1,name="",sequence = "",data={}):
        """Initialise a ChipInterval.

        Args:
            chr: Reference sequence name / chromosome.
            start: Start coordinate of the interval.
            end: End coordinate of the interval.
            strand: Strand indicator; defaults to '*' (unstranded).
            score: Probe or enrichment score; defaults to 0.0.
            readcount: Number of reads/probes; defaults to -1 (unset).
            name: Optional label for the interval; defaults to ''.
            sequence: Optional genomic sequence; defaults to ''.
            data: Optional dict of additional attributes; defaults to {}.
        """
        Interval.__init__(self, chr, start, end, strand=strand, score=score, readcount = readcount,name=name,sequence = sequence,data=data)
        self.parents = []
        self.children = []

    def addChild(self, child):
        """Add a child interval to this interval's children list.

        The child is only added if it is not already present.  A back-
        reference from the child to this interval is added to
        ``child.parents``.

        Args:
            child: A ChipInterval to add as a child of this interval.
        """
        #assert child not in self.children
        if child not in self.children:
            child.parents.append(self)
            self.children.append(child)

    def removeChild(self, child):
        """Remove a child interval from this interval's children list.

        Also removes the corresponding back-reference from ``child.parents``.
        The correctness of this method has not been fully verified.

        Args:
            child: The ChipInterval to remove from ``self.children``.
        """
        child.parents.remove(self)
        self.children.remove(child)

    def childScores(self):
        """Return the score attribute of each child interval.

        Returns:
            A list of score values, one per element in ``self.children``,
            in the same order as ``self.children``.
        """
        return [x.score for x in self.children]

    def childAvg(self):
        """Placeholder for computing the average score across child intervals.

        Not yet implemented.
        """
        pass

    def childMedian(self):
        """Placeholder for computing the median score across child intervals.

        Not yet implemented.
        """
        pass

    def makeValMap(self,value = 'readcount'):
        """Build a per-base value map by averaging child interval attributes.

        Creates ``self.valMap``, a numpy array of length ``len(self)``
        initialised to -1.  For each base position covered by at least one
        child interval the stored value is the mean of the specified
        attribute across all children that cover that base.

        Note: An alternative implementation exists in a commented-out block
        in the source; both approaches are noted as unverified.

        Args:
            value: Name of the attribute on each child ChipInterval whose
                values are averaged.  Defaults to ``'readcount'``.
        """
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
        """Plot probe scores across this interval using rpy2.

        Opens an X11 window and draws a step-style line plot.  Each child
        probe is drawn as a horizontal segment at its score level spanning
        its start to end coordinates.  If ``self.valMap`` has not yet been
        computed, ``makeValMap`` is called automatically.
        """
        if 'valMap' not in self.__dict__:
            self.makeValMap()
        robjects.r.x11()
        #robjects.r.plot(range(self.start,self.end+1),self.valMap,ylab="",type="l",lwd=2,main=str(self))
        robjects.r.plot((self.children[0].start,self.children[0].end),(self.children[0].score,self.children[0].score),type="l",lwd = 2,ylim=(min(c.score for c in self.children),max(c.score for c in self.children)))
        for x in self.children[1:]:
            robjects.r.lines((x.start,x.end),(x.score,x.score),lwd=2)

    def plot(self):
        """Convenience wrapper that calls plotVals to display the interval.

        Equivalent to calling ``self.plotVals()`` directly.
        """
        self.plotVals()

#    def uniqifySig(self):
#        keys = {}
#        for e in self.significant:
#            keys[e] = 1
#        self.significant = keys.keys()

    def scan(self,permuted,windowSize,threshold):
        """Scan child probes with a sliding window to identify significant regions.

        Sorts ``self.children`` in place and slides a window of
        ``windowSize`` probes across them.  For each window, computes the
        mean probe score and compares it against a pre-computed permutation
        distribution to obtain an empirical p-value.  Probes in windows
        whose p-value is at or below ``threshold`` are added to
        ``self.significant``.

        Args:
            permuted: A dict keyed by window size whose values are numpy
                arrays of maximum-window-mean values from permuted data (as
                produced by ``getRandomDist``).  The key ``windowSize`` must
                be present.
            windowSize: Number of consecutive probes in each sliding window.
            threshold: Maximum empirical p-value (proportion of permuted
                values >= observed mean) for a window to be considered
                significant.
        """
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
    """Container for one NimbleGen array's worth of probe data.

    Deprecated — this class is marked for deletion in the source.

    Parses a NimbleGen GFF file on construction and organises the resulting
    ChipInterval probe objects by chromosome.

    Attributes:
        fname: Path to the NimbleGen GFF file that was parsed.
        sampleName: Human-readable label for this sample.
        probeData: Dict mapping chromosome name to a list of ChipInterval
            objects for probes on that chromosome.
    """

    def __init__(self, fname, sampleName):
        """Initialise a ChipData container by parsing a NimbleGen GFF file.

        Args:
            fname: Path to the NimbleGen GFF output file to parse.
            sampleName: Label for this array sample.
        """
        self.fname = fname
        self.sampleName = sampleName
        self.probeData = {}

        #Populate self.probeData
        ChipIter = parseNimblegen(fname)
        for ci in ChipIter:
            if ci.chr not in list(self.probeData.keys()):
                self.probeData[ci.chr] = []
            self.probeData[ci.chr].append(ci)

    def sort(self):
        """Sort probe lists for all chromosomes in place.

        Iterates over ``self.data`` (note: the attribute populated on
        construction is ``self.probeData``; this method references
        ``self.data`` which may not exist).
        """
        for k in self.data.keys():
            self.data[k].sort()

    def shuffle(self,chr):
        """Shuffle probe scores for a chromosome in place.

        Note: This method is not yet correctly implemented — ``random.shuffle``
        operates on the temporary ``vals`` list and does not modify
        ``self.probeData``.

        Args:
            chr: Chromosome key to look up in ``self.probeData``.

        Returns:
            None (``random.shuffle`` always returns None).
        """
        vals = [x.score for x in self.probeData[chr]]
        return random.shuffle(vals)

#End crap

def nimblegenIter(fname):
    """Yield ChipInterval objects parsed from a NimbleGen GFF output file.

    Skips comment lines (starting with '#') and extracts chromosome,
    start, end, score, and probe name from each data row.

    Args:
        fname: Path to a NimbleGen GFF file.

    Yields:
        ChipInterval objects, one per non-comment line in the file.
    """
    handle = open(fname,'r')
    for line in handle:
        if line.startswith("#"): continue
        tokens = line.split("\t")
        pname = tokens[8].split(";")[1].split("=")[1]
        yield ChipInterval(tokens[0],tokens[3],tokens[4],score=tokens[5],name=pname)

def parseNimblegen(fname):
    """Parse an entire NimbleGen GFF file into a list of ChipInterval objects.

    Convenience wrapper around ``nimblegenIter`` that collects all intervals
    into a list rather than lazily yielding them.

    Args:
        fname: Path to a NimbleGen GFF file.

    Returns:
        A list of ChipInterval objects, one per non-comment line in the
        file.
    """
    iter = nimblegenIter(fname)
    rtrn = []
    for i in iter:
        rtrn.append(i)
    return rtrn

def joinNimblegenIntervals(intervals,start='start',end='end',offset=1000):
    """Merge overlapping NimbleGen probe intervals into enriched regions.

    Sorts the probe list and iterates through it, merging any probes that
    intersect (with optional extension by ``offset``) into a single
    ChipInterval.  Each merged interval stores its constituent probes as
    children and resets its name and score.

    Returns the input list unchanged if it is empty.

    Args:
        intervals: A list of ChipInterval objects (typically from
            ``parseNimblegen``).  The list is sorted in place.
        start: Attribute name used as the start coordinate when testing
            for intersection.  Defaults to ``'start'``.
        end: Attribute name used as the end coordinate when testing for
            intersection.  Defaults to ``'end'``.
        offset: Number of bases by which each interval is extended before
            testing for overlap, effectively merging probes within this
            distance.  Defaults to 1000.

    Returns:
        A list of merged ChipInterval objects representing independent
        enriched regions, each with a ``children`` list of the constituent
        probe intervals.
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
    """Extract scores from a list of probe intervals into a numpy array.

    Args:
        probes: A list of ChipInterval (or any object with a ``score``
            attribute) objects.

    Returns:
        A numpy array of dtype float32 containing the score of each probe
        in the same order as the input list.
    """
    return np.array([x.score for x in probes],dtype='f')

def getRandomDist(probes,nRandom,windowSize):
    """Build an empirical null distribution of maximum sliding-window means.

    Repeatedly shuffles the probe score array in place, slides a window of
    ``windowSize`` across it, records the maximum window mean for each
    shuffle, and returns all maxima as a numpy array.  This distribution is
    used to compute empirical p-values for observed window means.

    Args:
        probes: A numpy array (or list) of numeric probe scores.  The array
            is shuffled in place during this function — pass a copy if the
            original order must be preserved.
        nRandom: Number of shuffle iterations (i.e. length of the returned
            distribution array).
        windowSize: Number of consecutive probes in each sliding window.

    Returns:
        A numpy array of dtype float32 with ``nRandom`` elements, where each
        element is the maximum window mean observed in one shuffle iteration.
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
    """Count permuted values at least as extreme as the observed score.

    Note: This function is not yet correctly implemented.

    Args:
        segScores: The observed test statistic (scalar or array) to compare
            against the permuted distribution.
        permuted: A numpy array of values from the null distribution (e.g.
            as returned by ``getRandomDist``).
        windowSize: Window size used to generate the distribution (not
            currently used in the comparison).

    Returns:
        The number of elements in ``permuted`` that are >= ``segScores``.
    """
    return len(permuted[permuted>=segScores])


def main():
    """Run the default ChIP-chip analysis pipeline.

    Discovers all ``.gff`` files in the current directory, loads and
    normalises them via ``continuousData.SimpleChIPData``, merges adjacent
    probes, and generates permutation-based null distributions for a set of
    predefined window sizes (5, 7, 9, 11 probes).  The resulting
    distributions are stored in ``permuted`` keyed by window size.
    """
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
