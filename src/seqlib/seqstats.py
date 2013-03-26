#!/usr/bin/env python
import math,prob,misc,sys
import numpy
import mySam
import pysam
import intervallib
import scipy.stats
from RNASeq.misc import rstrips
import getopt
#from rpy2 import robjects
#from seqtools.genome import chr_lengths,genome_length

"""Collection of utilities for determining peak enrichment in xxx-Seq experiments"""

#################
#Main
#################
def main():
    expBam = sys.argv[1]
    ctlBam = sys.argv[2]
    bedFile = sys.argv[3]
    smRNApeakSeq(expBam,ctlBam,bedFile,filter=False,useStrand=useStrand)
    return


#########################
#Wrappers
########################
def smRNApeakSeq(expBam,ctlBam,bedFile,cutoff = 0.0001,filter=True,useStrand=True):
    #open files
    expHandle = pysam.Samfile(expBam,'rb')
    ctlHandle = pysam.Samfile(ctlBam,'rb')
    
    #Get normalization factor
    sys.stderr.write("Segmenting genome for Experimental BAM %s ...\n" % expBam)
    expBins = getSegmentCounts(expHandle)
    sys.stderr.write("Segmenting genome for Control BAM %s ...\n" % ctlBam)
    ctlBins = getSegmentCounts(ctlHandle)
    
    sys.stderr.write("Selecting non-zero indices ...\n")
    index = getNonZeroIndices(expBins,ctlBins)
    sys.stderr.write("Determining normalization factor ...\n")
    alpha = getAlpha(expBins,ctlBins,index)
    
    sys.stderr.write("alpha = %.4f\n" % alpha)
    
    del expBins
    del ctlBins
    del index
    
    #Loop over intervals
    sys.stderr.write("Testing intervals in %s...\n" % bedFile)
    results=[]
    bedIter = intervallib.parseBed(bedFile)
    for bed in bedIter:
        if useStrand:
            pVal,nExp,nCtl = testInterval(bed,expHandle,ctlHandle,alpha)
        else:
            pVal,nExp,nCtl = testIntervalNoStrand(bed,expHandle,ctlHandle,alpha)
        bed.data['pVal'] = pVal
        bed.data['nExp'] = nExp
        bed.data['nCtl'] = nCtl
        results.append(bed)
        
    #Correct for multiple tests
    #(Benjamini-Hochberg)
    sys.stderr.write("Correcting for multiple tests (%d)...\n" % len(results))
    results=multipleTestingCorrection(results)
    
    #Ran    k order by ascending q-value
    qVals = [x.data['qVal'] for x in results]
    qValRanks = misc.rank(qVals)
    
    sys.stderr.write("Printing results for %d tests..." % len(qValRanks))
    
    #Print header
    print "#chr\tstart\tend\tname\tscore\tstrand\tpVal\tqVal\tnExp\tnCtl"
    
    #This takes forever
    #count = 0
    #for i in xrange(len(qValRanks)):
    #    count += 1
    #    if count % 1000 == 0:
    #        sys.stderr.write("%g\n" % count)
    #    pos = qValRanks.index(i)
    #    res = results[pos]
    #    if not filter:
    #        print res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl'])
    #    else:
    #        if res.data['qVal'] <= cutoff:
    #           print res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl'])
    #sys.stderr.write("Done!\n")
    #return
    
    #Rank ordering output is too slow...just output and filter later.
    count = 0
    for res in results:
        count +=1
        if count % 1000 == 0:
            sys.stderr.write("%g\n" % count)
        if not filter:
            print res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl'])
        else:
            if res.data['qVal'] <= cutoff:
                print res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl'])
    sys.stderr.write("Done!\n")
    return
    
####################
#Normalization Functions
####################

def normDiff(expSum,ctlSum):
    """Takes the sum of reads (expSum) within a given interval and subtracts the sum of reads from the
    input or isotype control (ctlSum) for the same interval and then divides by the sqrt(expSum) to adjust for variance:
        (expSum-ctlSum)/sqrt(expSum)

    """ 
    return (expSum-ctlSum)/math.sqrt(expSum)

#####################
#Statistical Functions
#####################
'''
Deprecated in favor of using scipy and/or rpy2 (soooo much faster...)

def cumBinom(nExp,adjCtl,P=0.5):
    """Calculates P-value from the cumulative distribution function for the binomial distribution,
    which corresponds to summing the tail of the distribution.
    nCtl = number of reads in region i from control sample (e.g. input DNA or IgG)
    nExp = number of reads in region i from experimental sample
    k = adjCtl = nCtl*slope of linear regression # to normalize control sample to exp
    n = k + nExp
    P = 0.5 #probability under null hypothesis thaqt tags should occur with equal liklihood from sample and control
    """
    adjCtl = int(adjCtl)
    Pval = 0.0
    for j in range(0,adjCtl+1):
        Pval+=prob.binomial(P,j,nExp+adjCtl)
    return Pval
'''

def cumBinom(nExp,adjCtl,P=0.5):
    """
    The expected frequency of normalized reads for a given bin is p=0.5, therefore there is an equal likelihood that a read
    will be from either the experimental or control sample. This function uses scipy.stats.binom to return the probability 
    of observing >= nExp ( ie. 1-Pr(X <= x) ) reads from a given bin where k = nExp+adjCtl and P=0.5
    """
    return 1-scipy.stats.binom.cdf(nExp-1,nExp+adjCtl,P)

def testInterval(interval,expHandle,ctlHandle,alpha):
    """
    #TODO:Make sure that this is only grabbing the appropriate strand and not both....this can be dangerous
    """
    
    #expCounter = mySam.Counter()
    expCounter = mySam.StrandCounter()
    #ctlCounter = mySam.Counter()
    ctlCounter = mySam.StrandCounter()
    expFetch = expHandle.fetch(interval.chr,interval.start,interval.end,callback=expCounter)
    ctlFetch = ctlHandle.fetch(interval.chr,interval.start,interval.end,callback=ctlCounter)
    
    if interval.isPlus():
        nExp,nCtl = expCounter.plusCount,ctlCounter.plusCount

    elif interval.isMinus():
        nExp,nCtl = expCounter.minusCount,ctlCounter.minusCount
    #print nExp
    #print nCtl
    return cumBinom(nExp,nCtl*alpha),nExp,nCtl*alpha

def testIntervalNoStrand(interval,expHandle,ctlHandle,alpha):
    expCounter = mySam.Counter()
    ctlCounter = mySam.Counter()
    expFetch = expHandle.fetch(interval.chr,interval.start,interval.end,callback=expCounter)
    ctlFetch = ctlHandle.fetch(interval.chr,interval.start,interval.end,callback=ctlCounter)
    
    nExp,nCtl = expCounter.mCounts,ctlCounter.mCounts
    
    return cumBinom(nExp,nCtl*alpha),nExp,nCtl*alpha

def multipleTestingCorrection(testedIntervals):
    """Takes a list of Pvals and returns Qvals according to Benjamini-Hochberg:
    Count = number of hypothesis tests
    Q-value = P-value*(Count/Rank)
    """
    pVals = [x.data['pVal'] for x in testedIntervals]
    num = len(pVals)
    ranks = misc.rank(pVals)
    for i in range(0,num):
        qVal = pVals[i]*(float(num)/(ranks[i]+1))
        testedIntervals[i].data['qVal'] = qVal
    return testedIntervals

def getLambda(nReads,readLength,searchSize=3080419480):
    """A set of randomly located mapped DNA/RNA fragments is equivalent to a global coverage level lambda, 
    whose value is the product of the number and mean length of mapped fragments divided by the mappable 
    search space length (genome size).
    
    returns lambda: a measure of expected coverage per base of the search space
    """
    
    return (nReads*readLength)/(float(searchSize))

def poissonProb(lamb,height):
    """
    ***THIS IS WRONG***
    I think that the correct lambda should be the per-base expectancy * the size of the peak, but I will have to check
    
    TODO:Currently does naive calculation of cdf by summing point probabilities (will fix that)
    
    Given a lambda value, the probability of observing a peak with a height >= H 
    is given by a sum of Poisson probabilities (1-cdf(height-1,lambda))
    
    Returns 1-cumulative density function = probability of finding a peak of height
    H or greater given a global per-base coverage value of k (assuming random background)
    """
    probs = 0.0
    for k in range(0,height-1):
        probs += ((math.e**(-lamb)*lamb**k)/prob.factorial(k))
    
    return 1-probs
    
    """
    OR
    return scipy.stats.poisson.cdf(height-1,lamb)
    
    """    
    

#########################
#Normalization utilities
#########################

def slope(xarray,yarray):
    """Uses numpy, in fact assumes that the list arguments are numpy arrays."""
    n = float(len(xarray))
    m = (n*sum(xarray*yarray)-sum(xarray)*sum(yarray))/(n*sum(xarray**2)-(sum(xarray))**2)
    return m

def intercept(xarray,yarray):
    """Uses numpy, in fact assumes that the list arguments are numpy arrays."""
    m = slope(xarray,yarray)
    n = float(len(xarray))
    b = (sum(yarray)-m*(sum(xarray)))/n
    return b

def getSegmentCounts(bamHandle,segSize=10000):
    chrs = bamHandle.references
    chr_lengths = bamHandle.lengths
    bins = numpy.zeros(sum(chr_lengths)/segSize+len(chrs))
    index = 0
    for x in xrange(0,len(chrs)):
        sys.stderr.write(chrs[x]+"\n")
        for i in xrange(0,chr_lengths[x],segSize):
            c = mySam.Counter()
            bamHandle.fetch(chrs[x],i,i+segSize,callback=c)
            bins[index] += (c.mCounts)
            index+=1
            #print c.mCounts
    return bins

def getNonZeroIndices(bins1,bins2):
    set1 = set(numpy.nonzero(bins1)[0])
    set2 = set(numpy.nonzero(bins2)[0])
    return list(set1.intersection(set2))

def getAlpha(expBins,ctlBins,index):
    return slope(ctlBins[index],expBins[index])

def getAlphaFromLinReg(exp,ctl,r):
    """Takes two lists of counts for paired regions (experimental and control) and performs a linear regression to normalize.
    Returned value is slope of linear regression to be used as scaling factor for control sample
    ****Deprecated****
    """

    r.plot(ctl,exp,ylab = "Experimental Counts", xlab = "Control Counts", pch=20)
    lm = r.lsfit(ctl,exp)
    gradient = lm['coefficients']['X']
    return gradient

########
#Main
########
help_message = '''
This is my own implementation of the Peak-Seq peak calling application

Usage:
python seqstats.py [options] -e <experimental.bam> -c <control.bam> -b <experimental.bed> >outfile.out

Options:
    -e | --expBam        Sorted and indexed bam file of the experimental sample
    -c | --ctlBam        Sorted and indexed bam file of the control sample
    -b | --expBed        Bed file of contiguous intervals from --expBam
    -s | --ignoreStrand  Ignore strand information when counting reads from each interval
    -h | --help          This helpful help message
    -v | --verbose       Verbose 
    -o | --outFile       Where to write the output
    --cutoff             Q-value cutoff (default: 0.0001)
    --filter             Filter output to only show results with Q-value greater than cutoff (default: off)
    
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def newMain(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(argv[1:], "he:c:b:o:sftv", ["help", "expBam=","ctlBam=","expBed=","output=","ignoreStrand","filter","cutoff","verbose="])
        except getopt.error, msg:
            raise Usage(msg)
        #Defaults
        verbose = False
        useStrand = True
        filter = False
        cutoff = 0.0001
        #option processing
        for option, value in opts:
            if option in ("-v","--verbose"):
                verbose = True
            if option in ("-e","--expBam"):
                expBam = value
            if option in ("-c","--ctlBam"):
                ctlBam = value
            if option in ("-b","--expBed"):
                expBed = value
            if option in ("-h","--help"):
                raise Usage(help_message)
            if option in ("-o","--output"):
                #output doesn't work quite yet
                outFile = value
            if option in ("-s","--ignoreStrand"):
                useStrand = False
            if option == "--cutoff":
                cutoff = float(value)
            if option == "--filter":
                filter = True

#        if outFile == None:
#            outFile = rstrips(expBed,".bed")+".out"    
        #Call Main with arguments
        smRNApeakSeq(expBam,ctlBam,expBed,filter=filter,cutoff=cutoff,useStrand=useStrand)
    except Usage,err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    return

if __name__ == "__main__":
    newMain()