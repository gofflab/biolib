#!/usr/bin/env python
"""Statistical utilities for peak enrichment analysis in RNA immunoprecipitation and ChIP-Seq experiments.

Implements a PeakSeq-like approach for comparing experimental (RIP or ChIP)
BAM files against input/IgG control BAM files.  The pipeline:

1. Segments the genome into fixed-size bins and counts reads in each bin for
   both the experimental and control samples.
2. Determines a global normalisation factor (alpha) via linear regression on
   bins that have reads in both samples.
3. Tests each interval in a BED file using a binomial model (reads from the
   experimental sample vs. alpha-scaled control reads) to assign p-values.
4. Corrects for multiple testing using Benjamini-Hochberg FDR correction.
5. Outputs results to stdout sorted or filtered by q-value.
"""
import getopt
import math
import sys

import numpy
import pysam
import scipy.stats

from . import intervallib, misc, mySam, prob

#from rpy2 import robjects
#from seqtools.genome import chr_lengths,genome_length

#################
#Main
#################
def main():
    """Legacy command-line entry point — reads three positional arguments and runs smRNApeakSeq.

    Expects sys.argv to contain: expBam ctlBam bedFile.  Calls smRNApeakSeq
    with filter=False and the module-level useStrand variable.  Prefer
    newMain() for proper option parsing.
    """
    expBam = sys.argv[1]
    ctlBam = sys.argv[2]
    bedFile = sys.argv[3]
    smRNApeakSeq(expBam,ctlBam,bedFile,filter=False,useStrand=useStrand)
    return


#########################
#Wrappers
########################
def smRNApeakSeq(expBam,ctlBam,bedFile,cutoff = 0.0001,filter=True,useStrand=True):
    """Runs the full smRNA/RIP-Seq peak-calling pipeline and writes results to stdout.

    Segments the genome, computes a normalisation factor between experimental
    and control BAM files, tests each BED interval with a binomial model,
    applies Benjamini-Hochberg FDR correction, and prints tab-delimited
    output.

    Args:
        expBam: Path to a sorted, indexed BAM file from the experimental
            (RIP/ChIP) sample.
        ctlBam: Path to a sorted, indexed BAM file from the control (IgG or
            input) sample.
        bedFile: Path to a BED file of candidate intervals to test.
        cutoff: Q-value threshold below which results are printed when filter
            is True (default: 0.0001).
        filter: If True, only print intervals with q-value <= cutoff.  If
            False, print all intervals (default: True).
        useStrand: If True, count only reads on the same strand as each
            interval.  If False, count all reads regardless of strand
            (default: True).
    """
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

    #Rank order by ascending q-value
    qVals = [x.data['qVal'] for x in results]
    qValRanks = misc.rank(qVals)

    sys.stderr.write("Printing results for %d tests..." % len(qValRanks))

    #Print header
    print("#chr\tstart\tend\tname\tscore\tstrand\tpVal\tqVal\tnExp\tnCtl")

    #This takes forever
    #count = 0
    #for i in range(len(qValRanks)):
    #    count += 1
    #    if count % 1000 == 0:
    #        sys.stderr.write("%g\n" % count)
    #    pos = qValRanks.index(i)
    #    res = results[pos]
    #    if not filter:
    #        print(res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl']))
    #    else:
    #        if res.data['qVal'] <= cutoff:
    #           print(res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl']))
    #sys.stderr.write("Done!\n")
    #return

    #Rank ordering output is too slow...just output and filter later.
    count = 0
    for res in results:
        count +=1
        if count % 1000 == 0:
            sys.stderr.write("%g\n" % count)
        if not filter:
            print(res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl']))
        else:
            if res.data['qVal'] <= cutoff:
                print(res.toBed()+"\t%g\t%g\t%d\t%d" % (res.data['pVal'],res.data['qVal'],res.data['nExp'],res.data['nCtl']))
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
    """Tests a single genomic interval for strand-aware read enrichment.

    Counts reads on the same strand as the interval from both the experimental
    and control BAM files, scales the control count by alpha, and returns
    a binomial p-value.

    Args:
        interval: An intervallib.Interval object with chr, start, end, and
            strand attributes.
        expHandle: A pysam AlignmentFile for the experimental sample.
        ctlHandle: A pysam AlignmentFile for the control sample.
        alpha: Normalisation factor (slope from getAlpha) used to scale
            control counts to match the experimental library size.

    Returns:
        A tuple (pVal, nExp, adjCtl) where pVal is the binomial p-value,
        nExp is the raw experimental read count, and adjCtl is the
        alpha-scaled control read count.
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
    """Tests a single genomic interval for read enrichment ignoring strand.

    Counts all reads (both strands) overlapping the interval from experimental
    and control BAM files, scales control count by alpha, and returns a
    binomial p-value.

    Args:
        interval: An intervallib.Interval object with chr, start, and end
            attributes.
        expHandle: A pysam AlignmentFile for the experimental sample.
        ctlHandle: A pysam AlignmentFile for the control sample.
        alpha: Normalisation factor used to scale control counts.

    Returns:
        A tuple (pVal, nExp, adjCtl) where pVal is the binomial p-value,
        nExp is the raw experimental read count, and adjCtl is the
        alpha-scaled control read count.
    """
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
    """Computes the slope of the ordinary least-squares regression line.

    Uses numpy arrays for efficient computation.  The slope is:
        m = (n*sum(x*y) - sum(x)*sum(y)) / (n*sum(x^2) - (sum(x))^2)

    Args:
        xarray: A numpy array of x (independent variable) values.
        yarray: A numpy array of y (dependent variable) values of the same
            length as xarray.

    Returns:
        The slope of the linear regression line (float).
    """
    n = float(len(xarray))
    m = (n*sum(xarray*yarray)-sum(xarray)*sum(yarray))/(n*sum(xarray**2)-(sum(xarray))**2)
    return m

def intercept(xarray,yarray):
    """Computes the y-intercept of the ordinary least-squares regression line.

    Uses numpy arrays for efficient computation.  The intercept is:
        b = (sum(y) - m*sum(x)) / n

    Args:
        xarray: A numpy array of x (independent variable) values.
        yarray: A numpy array of y (dependent variable) values of the same
            length as xarray.

    Returns:
        The y-intercept of the linear regression line (float).
    """
    m = slope(xarray,yarray)
    n = float(len(xarray))
    b = (sum(yarray)-m*(sum(xarray)))/n
    return b

def getSegmentCounts(bamHandle,segSize=10000):
    """Counts reads in fixed-size genomic bins across all chromosomes in a BAM file.

    Iterates over all chromosomes and divides each into bins of segSize base
    pairs, counting the total number of reads per bin using mySam.Counter.

    Args:
        bamHandle: A pysam AlignmentFile opened for reading.
        segSize: Bin size in base pairs (default: 10000).

    Returns:
        A numpy array of read counts, one element per bin, ordered by
        chromosome then genomic position.
    """
    chrs = bamHandle.references
    chr_lengths = bamHandle.lengths
    bins = numpy.zeros(sum(chr_lengths)//segSize+len(chrs))
    index = 0
    for x in range(0,len(chrs)):
        sys.stderr.write(chrs[x]+"\n")
        for i in range(0,chr_lengths[x],segSize):
            c = mySam.Counter()
            bamHandle.fetch(chrs[x],i,i+segSize,callback=c)
            bins[index] += (c.mCounts)
            index+=1
            #print c.mCounts
    return bins

def getNonZeroIndices(bins1,bins2):
    """Returns the indices of bins that have non-zero counts in both arrays.

    Used to restrict linear regression normalisation to bins that are
    informative in both the experimental and control samples.

    Args:
        bins1: A numpy array of read counts (e.g. experimental sample bins).
        bins2: A numpy array of read counts (e.g. control sample bins) of
            the same length as bins1.

    Returns:
        A list of integer indices where both bins1 and bins2 have non-zero
        values.
    """
    set1 = set(numpy.nonzero(bins1)[0])
    set2 = set(numpy.nonzero(bins2)[0])
    return list(set1.intersection(set2))

def getAlpha(expBins,ctlBins,index):
    """Computes the normalisation factor (alpha) between experimental and control samples.

    Fits a linear regression through the origin on the subset of bins
    specified by index, treating control counts as x and experimental counts
    as y.  The slope is used to scale the control sample to the experimental
    library size.

    Args:
        expBins: Numpy array of per-bin read counts for the experimental
            sample.
        ctlBins: Numpy array of per-bin read counts for the control sample.
        index: List of integer indices identifying informative bins (non-zero
            in both arrays).

    Returns:
        Alpha (float): the slope of the linear regression, used as the
        multiplicative scaling factor for control counts.
    """
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
    """Exception raised for command-line usage errors in seqstats.

    Attributes:
        msg: Human-readable explanation of the error or the help message.
    """
    def __init__(self, msg):
        """Initialises a Usage exception.

        Args:
            msg: Human-readable error or help text.
        """
        self.msg = msg

def newMain(argv=None):
    """Command-line entry point for the seqstats peak-calling pipeline.

    Parses command-line options and delegates to smRNApeakSeq.  Supports
    optional strand-specific counting, q-value filtering, and verbose output.

    Args:
        argv: List of command-line argument strings.  Defaults to sys.argv
            when None.

    Returns:
        2 on usage error, None on success.

    Raises:
        SystemExit: Indirectly via sys.exit() on usage error.
    """
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt(argv[1:], "he:c:b:o:sftv", ["help", "expBam=","ctlBam=","expBed=","output=","ignoreStrand","filter","cutoff","verbose="])
        except getopt.error as msg:
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
    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg), file=sys.stderr)
        print("\t for help use --help", file=sys.stderr)
        return 2
    return

if __name__ == "__main__":
    newMain()
