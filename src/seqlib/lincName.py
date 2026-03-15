#!/usr/bin/env python
"""Assigns systematic names to lincRNA loci based on proximity to RefSeq genes.

Implements the naming scheme described in Guttman et al. for long intergenic
non-coding RNA (lincRNA) loci:

- If the 5' end of a lincRNA overlaps the 5' end of a protein-coding gene on
  the opposite strand by less than the overlap threshold, the lincRNA is named
  'linc-<GENE>-BP' (bidirectional promoter).
- If a lincRNA overlaps any protein-coding gene on the opposite strand without
  satisfying the bidirectional criterion, it is named 'linc-<GENE>-AS'
  (antisense).
- Otherwise, the lincRNA is named after the nearest downstream protein-coding
  gene on the same strand: 'linc-<GENE>' (single lincRNA) or
  'linc-<GENE>-<N>' (multiple lincRNAs near the same gene).

Requires GTFlib, dbConn, and intervallib packages, and a connection to the
UCSC genome browser MySQL server.
"""

############
#Imports
############
import bisect
import copy
import getopt
import sys

import dbConn
import GTFlib
from misc import rstrips

############
#Constants
############
overlapThreshold = 0.20
extensionLength = 500 #grow 5'end of lincRNA by this many bases to test for Bidirectional promoter
strandLookup = {'+':'-','-':'+'}

help_message = '''
Created on Aug 27, 2010
@author: lgoff

Usage: python lincName.py [options] <gtfFile.gtf>

Options:
    -g | --genome  [Default : hg19]   Determines what build of the genome is used to fetch RefSeq transcripts
                    around which lincNames are chosen.
                    
    -h | --help       Displays this helpful help screen
    
    -v                Verbose
    
    -o | --output    [Default : <gtfFile_named.gtf>] Determines output file
'''

############
#Classes
############
class Usage(Exception):
    """Exception raised for command-line usage errors in lincName.

    Attributes:
        msg: Human-readable explanation of the error or the help message.
    """
    def __init__(self, msg):
        """Initialises a Usage exception.

        Args:
            msg: Human-readable error or help text.
        """
        self.msg = msg


############
#Functions
############

def test5PrimeOverlap(lincInt,geneInt):
    """Determines whether the overlap between a lincRNA and a gene is at the lincRNA 5' end.

    Tests whether a lincRNA interval overlaps a protein-coding gene such that
    the overlap is at the 5' end of the lincRNA (and also involves the 5' end
    of the gene on the opposite strand).  Used to identify bidirectional
    promoter pairs.

    Note: may not give correct results when a lincRNA completely spans a
    protein-coding gene on the opposite strand.

    Args:
        lincInt: An interval object for the lincRNA with strand, start, and
            end attributes.
        geneInt: An interval object for the overlapping protein-coding gene
            with strand, start, and end attributes.

    Returns:
        True if the overlap is at the 5' end of lincInt; False otherwise.

    Raises:
        AssertionError: If the two intervals do not overlap.
        ValueError: If the strand of lincInt cannot be determined.
    """
    assert lincInt.overlaps(geneInt)
    if lincInt.strand == "+":
        if lincInt.start <= geneInt.end and lincInt.end > geneInt.end:
            return True
        else:
            return False
    elif lincInt.strand == "-":
        if geneInt.start <= lincInt.end and geneInt.end > lincInt.end:
            return True
        else:
            return False
    else:
        raise ValueError("Could not determine")

def bpOverlap(lincInt,geneInt):
    """Returns the number of base pairs of overlap between two genomic intervals.

    Sorts the four boundary coordinates and computes the inner distance as the
    length of the shared region.

    Args:
        lincInt: An interval object with start and end attributes.
        geneInt: An interval object with start and end attributes that must
            overlap with lincInt.

    Returns:
        Integer number of overlapping base pairs between the two intervals.

    Raises:
        AssertionError: If the two intervals do not overlap.
    """
    assert lincInt.overlaps(geneInt), "%s and %s do not overlap" % (lincInt.name,geneInt.name)
    bounds = [lincInt.start,lincInt.end,geneInt.start,geneInt.end]
    bounds.sort()
    #range = bounds[3]-bounds[0]
    overlap = bounds[2]-bounds[1]
    return overlap

def printLincs(handle,lincs):
    """Writes a collection of lincRNA GTF records to a file handle.

    Args:
        handle: Writable file-like object to receive the GTF output.
        lincs: Iterable of lincRNA objects, each exposing a getGTF() method
            that returns a GTF-formatted string.
    """
    for linc in lincs:
        print(linc.getGTF(), end=' ', file=handle)

############
#Main
############

def main(gtfFile,genome='hg19'):
    """Assigns systematic names to all lincRNA loci in a GTF file.

    Reads lincRNA transcript models from gtfFile, retrieves protein-coding
    RefSeq transcripts for the specified genome build, and applies the
    bidirectional promoter, antisense, and proximity naming rules to produce
    a set of named lincRNA objects.

    Args:
        gtfFile: Path to a GTF file of unannotated lincRNA loci (as produced
            by Cufflinks or similar assemblers).
        genome: UCSC genome build identifier used to fetch RefSeq transcripts
            (default: 'hg19').

    Returns:
        A set of lincRNA gene objects with updated name attributes following
        the systematic naming convention.
    """
    #Parse GTF File for lincs
    lincIter = GTFlib.GTFGeneIterator(gtfFile,verbose=verbose)

    #Retrieve and index RefSeq genes
    refSeqs = dbConn.fetchRefSeqIntervalsIndexed(genome=genome,proteinCodingOnly=True,verbose=verbose)

    #Results container
    res = set([])

    #Container for gene:linc assoc.
    geneLincs = {}

    #Loop through lincRNAs
    for linc in lincIter:
        flag = False
        bdFlag = False #True if linc is bidirectional
        asFlag = False #True if linc is antisense
        #Convert to Interval
        interval = linc.toInterval()

        #Test for weird chromosome (ie. not in refSeqs.keys() )
        if interval.chr not in refSeqs.keys():
            res.add(linc)
            continue

        #Bug tracking only
        if verbose:
            sys.stderr.write(str(interval)+"\n")

        #Get list of gene positions that are relevant
        senseGeneStarts = [x.start for x in refSeqs[interval.chr][interval.strand]]
        senseGeneEnds = [x.end for x in refSeqs[interval.chr][interval.strand]]

        #Get opposite strand to test
        testStrand = strandLookup[interval.strand]

        #Test overlap with genes on opposite strand
        for gene in refSeqs[interval.chr][testStrand]:
            extendedInterval = copy.copy(interval)
            extendedInterval.grow5_prime(extensionLength)

            if extendedInterval.overlaps(gene):
                #If 5' end of linc overlaps the 5' of a coding gene on the opposite strand,
                #by more than 0bp but less than min(BP_THRESH * length(L), BP_THRESH * length(coding gene))
                #THEN name linc "linc-[HUGO_GENE_NAME]-BP"
                overlap = bpOverlap(extendedInterval,gene)
                fivePrime = test5PrimeOverlap(extendedInterval,gene)
                cutoff = min(len(extendedInterval)*overlapThreshold,gene.intervalLen()*overlapThreshold)
                if fivePrime and overlap <= cutoff:
                    linc.propogateLincName("linc-%s-BP" % gene.name)
                    linc.addAttribute("bidirectional_prom",gene.name)
                    res.add(linc)
                    flag = True
                    bdFlag = True
                    #break
                    continue

                #TODO FIX this so that ANY overlap that is not a BP becomes and -AS
                if not bdFlag:
                    linc.propogateLincName("linc-%s-AS" % gene.name)
                linc.addAttribute("antisense",gene.name)
                res.add(linc)
                flag = True
                asFlag = True
                break
        #ELSE find the closest coding gene on the same strand as the L, starting from the 3' end of the linc.
        #Suppose its HUGO name is NCG1.Add L to a list of lincs to be named after NCG1.
        if not flag:
            if interval.strand == "+":
                nearestGeneIdx = bisect.bisect(senseGeneStarts,interval.end) #choose most adjacent gene 3' to lincRNA
            elif interval.strand == "-":
                nearestGeneIdx = bisect.bisect(senseGeneEnds,interval.start)-1
            try:
                nearestGene = refSeqs[interval.chr][interval.strand][nearestGeneIdx]
            except IndexError:
                #If I cannot find the nearestGene (e.g. end of chromosome or something, just push linc to results
                #and deal with them later. (for now)

                #print nearestGeneIdx
                #print interval.toBed()
                res.add(linc)
                continue
            geneLincs.setdefault(nearestGene.name,[]).append(linc)

    #Evaluate container for linc:gene assocs
    """
    FOREACH coding gene G in the table above:
    IF there's only one linc to be named after G THEN
        name that linc "linc-G"
    ELSE
        sort the list of lincs by proximity to G, with the closest linc at the front of the list
        FOR i = 1 to #number of lincs named after G
            name linc i "linc-G-i"
    """
    for k,v in geneLincs.iteritems():
        if len(v) == 1:
            v[0].propogateLincName("linc-%s" % (k))
            res.add(v[0])
        elif len(v) >1:
            if v[0].strand == "+":
                v.sort(reverse=True)
            elif v[0].strand == "-":
                v.sort()
            for i in xrange(len(v)):
                v[i].propogateLincName("linc-%s-%d" % (k,i+1))
                res.add(v[i])
    return res

############
#Tests
############
def test():
    """Runs a full naming test using hardcoded Broad Institute file paths.

    Calls main() on a hard-coded lincRNA GTF file, writes the named output
    to a companion file, and prints a completion message to stderr.  Intended
    for interactive/development testing only.
    """
    fname = '/seq/rinnscratch/cole/ftp/assemblies/linc_catalog.gtf'
    outHandle = open('/seq/rinnscratch/cole/ftp/assemblies/linc_catalog_named.gtf','w')
    verbose=True
    lincs = main(fname)
    printLincs(outHandle,lincs)
    sys.stderr.write("Done!"+"\n")
    return



############
#Orders
############
if __name__=="__main__":
    #test()
    argv = sys.argv
    #default settings
    genome = "hg19"
    verbose = False
    outFile = None
    try:
        try:
            opts,args = getopt.getopt(argv[1:],"hg:o:v",["help","genome","output"])
        except getopt.error as msg:
            raise Usage(msg)

        #option processing
        for option,value in opts:
            if option in ("-g","--genome"):
                genome = value
            if option in ("-h","--help"):
                raise Usage(help_message)
            if option == "-v":
                verbose = True
            if option in ("-o","--output"):
                outFile = value

        #debugging
        #print opts
        #print args

        try:
            assert len(args)==1
            gtfFile = args[0]
        except:
            raise Usage(help_message)
        baseName = rstrips(gtfFile,".gtf")
        if verbose:
            sys.stderr.write("Naming lincs in file %s using RefSeq transcripts in genome %s.\n" % (gtfFile,genome))
        lincs = main(gtfFile,genome=genome)
        if outFile == None:
            outFile = (baseName+"_named.gtf")
        if verbose:
            sys.stderr.write("Writing output to %s.\n" % outFile)
        outHandle = open(outFile,'w')
        printLincs(outHandle,lincs)
        if verbose:
            sys.stderr.write("Done!\n")
    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg), file=sys.stderr)
        sys.exit()

