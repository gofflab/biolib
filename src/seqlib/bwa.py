'''
Python wrappers for the BWA short-read alignment algorithm.

Provides helper functions for submitting BWA alignment jobs to an LSF
cluster (``bsub``), converting SAM output to sorted BAM files, and parsing
SAM records.  Also includes utilities for converting pileup output to UCSC
wiggle format.

The module-level ``prefix`` and ``ref_index`` constants point to the hg18
reference genome used by the original author; update these for other
references.

Created on Jul 30, 2009

@author: lgoff

Example of base commands:
BWA Align:
    bwa aln -c /seq/compbio-hp/lgoff/genomes/hg18/hg18.fa test.fastq >test.sai
BWA SAMSE:
     bwa samse /seq/compbio-hp/lgoff/genomes/hg18/hg18.fa test.sai test.fastq
'''
import copy
import os

from .Alignment import *

prefix = "/seq/compbio-hp/lgoff/genomes/hg18/hg18.fa"
ref_index = prefix+".fai"

#=================
class SAMAlignment(Alignment):
    """SAM alignment record with CIGAR and quality-string fields.

    Extends the Alignment base class with the two SAM-specific fields that
    are not part of the generic Alignment interface.

    Attributes:
        qual: ASCII-encoded base-quality string (SAM field 11).
        cigar: CIGAR string describing the alignment operations (SAM field 6).
    """

    def __init__(self,readname,chr,start,end,strand,score,readcount,readsequence,cigar,qualstring):
        """Initialise a SAMAlignment.

        Args:
            readname: Query template name (SAM field 1).
            chr: Reference sequence name / chromosome (SAM field 3).
            start: 1-based leftmost mapping position (SAM field 4).
            end: Computed end position (start + read length - 1).
            strand: Strand of the alignment, '+' or '-'.
            score: Mapping quality score (SAM field 5).
            readcount: Number of reads represented (typically 1).
            readsequence: Read sequence bases (SAM field 10).
            cigar: CIGAR string (SAM field 6).
            qualstring: ASCII-encoded base-quality string (SAM field 11).
        """
        Alignment.__init__(self,readname,chr,start,end,strand,score=readcount,readcount = readcount,readsequence=readsequence)
        self.qual = qualstring
        self.cigar = cigar

def SAMReader(fname):
    """Iterate over SAM alignment records from a file.

    Args:
        fname: Path to the SAM file.

    Yields:
        An Interval object for each alignment record in the file.
    """
    handle = open(fname,'r')
    for line in handle:
        aln = parseSAMString(line)
        yield aln.toInterval()

def parseSAMString(samstring):
    """Parse a single SAM-format line into a SAMAlignment object.

    The end position is derived from the start position plus the length of
    the read sequence field; this is only correct for non-spliced alignments.

    Args:
        samstring: A single tab-delimited SAM record line (trailing whitespace
            is stripped internally).

    Returns:
        A SAMAlignment instance populated from the SAM fields.
    """
    tokens = samstring.rstrip().split("\t")
    readname = tokens[0]
    chr = tokens[2]
    start = int(tokens[3])
    end = int(tokens[3])+len(tokens[9])-1
    strand = strandFlag(tokens[1])
    score = int(tokens[4])
    readcount = 1
    readsequence = tokens[9]
    cigar = tokens[5]
    qualstring = tokens[10]
    return SAMAlignment(readname,chr,start,end,strand,score,readcount,readsequence,cigar,qualstring)

def joinSAMIntervals(iter,start='start',end='end',offset=0):
    """Merge overlapping SAM intervals into non-overlapping intervals, per strand.

    Groups intervals by strand ('+' or '-'), then iterates through each
    group in order and merges any pair of intervals that intersect (with
    optional extension by ``offset``).  Each merged interval stores its
    constituent child intervals and reports their count as ``readcount``.

    The SAM file must be sorted with ``samtools sort`` before use.

    Args:
        iter: An iterable of Interval (or Alignment) objects already
            loaded from a sorted SAM file.  Each must have a ``strand``
            attribute of '+' or '-'.
        start: Name of the start-coordinate attribute used when testing
            intersection.  Defaults to 'start'.
        end: Name of the end-coordinate attribute used when testing
            intersection.  Defaults to 'end'.
        offset: Number of bases by which interval extents are extended
            before testing for overlap.  Defaults to 0.

    Returns:
        A dict with keys '+' and '-', each mapping to a list of merged
        Interval objects for that strand.  Each merged interval has a
        ``readcount`` equal to the number of constituent child reads and
        a ``children`` list of those child intervals.
    """

    overlapping_plus = []
    overlapping_minus = []
    for interval in iter:
        if interval.strand=="+":
            overlapping_plus.append(interval)
        elif interval.strand=="-":
            overlapping_minus.append(interval)
        else:
            continue
    res = {}
    for i in ("+","-"):
        print(i)
        if i =="+":
            intervals = overlapping_plus
        elif i =="-":
            intervals = overlapping_minus
        else: continue # should not have to resort to this
        non_overlapping = []
        current = copy.copy(intervals[0])
        current.addChild(copy.copy(current))
        current.readcount = -1
        for x in intervals[1:]:
            next = copy.copy(x)
            if current.intersects(next,start=start,end=end,offset=offset):
                current.end = max(current.end,next.end)
                current.addChild(copy.copy(next))
            else:
                current.readcount=len(current.children)
                non_overlapping.append(current)
                #print current
                current = next
                current.addChild(copy.copy(current))
                current.readcount=-1
        current.readcount = len(current.children)
        non_overlapping.append(current)
        res[i] = non_overlapping
    return res

def bwaAlignSubmit(files,mismatches=2,queue='hugemem'):
    """Submit BWA alignment jobs (``bwa aln``) to an LSF cluster.

    For each input FASTQ file, constructs and submits an LSF ``bsub`` job
    that runs ``bwa aln`` against the module-level ``prefix`` reference and
    writes a ``.sai`` alignment index file.

    Args:
        files: A list of FASTQ file paths to align.
        mismatches: Maximum number of mismatches allowed in the seed region
            (passed to ``bwa aln -n``).  Defaults to 2.
        queue: LSF queue name to submit jobs to.  Defaults to 'hugemem'.
    """
    for fname in files:
        shortname = fname.rstrip(".fastq")
        command = "bsub -q %s -N -o /dev/null -P BWA_Align 'bwa aln -c -n %d %s %s >%s.sai 2>%s.e'" % (queue,mismatches,prefix,fname,shortname,shortname)
        os.system(command)
    return

def bwaSamseSubmit(files,mismatches=2,queue='broad'):
    """Submit BWA SAM conversion jobs (``bwa samse``) to an LSF cluster.

    For each ``.sai`` file, constructs and submits an LSF ``bsub`` job that
    runs ``bwa samse`` to convert the alignment index back to SAM format,
    writing a ``.sam`` file.  Assumes a matching ``.fastq`` file exists with
    the same base name.

    Args:
        files: A list of ``.sai`` file paths to convert.
        mismatches: Unused parameter kept for interface compatibility.
            Defaults to 2.
        queue: LSF queue name to submit jobs to.  Defaults to 'broad'.
    """
    for fname in files:
        shortname = fname.rstrip(".sai")
        command = "bsub -q %s -N -o /dev/null -P BWA_Samse 'bwa samse %s %s.sai %s.fastq >%s.sam 2>%s.e'" % (queue,prefix,shortname,shortname,shortname,shortname)
        os.system(command)
    return

def makeBam(files,queue='broad'):
    """Submit SAM-to-BAM conversion jobs (``samtools view``) to an LSF cluster.

    For each SAM file, constructs and submits an LSF ``bsub`` job that uses
    ``samtools view`` to convert it to a BAM file indexed against the
    module-level ``ref_index`` FASTA index.

    Args:
        files: A list of SAM file paths to convert.
        queue: LSF queue name to submit jobs to.  Defaults to 'broad'.
    """
    for fname in files:
        shortname = fname.rstrip("*.sam")
        command = "bsub -q %s -N -o /dev/null -P SAM2BAM 'samtools view -h -bt %s -o %s.bam %s 2>%s.bam.e'" % (queue,ref_index,shortname,fname,shortname)
        os.system(command)
    return

def samSort(files,queue='broad'):
    """Sort BAM files by coordinate using ``samtools sort``.

    Iterates over a list of BAM files, printing a status message for each,
    and runs ``samtools sort`` locally (not via LSF) to produce a
    ``*_sorted.bam`` output file.

    Args:
        files: A list of BAM file paths to sort.
        queue: Unused parameter kept for interface consistency with other
            submit functions.  Defaults to 'broad'.
    """
    for fname in files:
        shortname = fname.rstrip("*.bam")+"_sorted"
        command = "samtools sort %s %s" % (fname,shortname)
        print("Sorting file: %s" % fname)
        os.system(command)
    return



def pileup2wig(fname,shortname,outDir=os.getcwd()+"/"):
    """Convert a samtools pileup file to strand-specific wiggle files.

    Reads a samtools pileup output file and writes two variableStep wiggle
    files: one for the plus strand (forward reads, '.' characters) and one
    for the minus strand (reverse reads, ',' characters).

    Args:
        fname: Path to the samtools pileup file to read.
        shortname: Base name used for the wiggle track labels and the output
            file names (``<shortname>_plus.wig`` and
            ``<shortname>_minus.wig``).
        outDir: Directory in which the output wiggle files are written.
            Defaults to the current working directory.
    """
    handle = open(fname,'r')
    preRef = ''
    prePos = -1
    prePlus = 0
    preMinus = 0

    plusHand = open(outDir+shortname+"_plus.wig",'w')
    minusHand = open(outDir+shortname+"_minus.wig",'w')

    def wigHeader(shortname,strand):
        """Build a UCSC wiggle track-definition header line.

        Args:
            shortname: Base name used in the track name and description.
            strand: Strand of the track; '+' produces a blue track,
                '-' produces a red track.

        Returns:
            A wiggle track header string suitable for the first line of a
            wiggle file.
        """
        if strand=="+":
            color = '0,0,255'
            sName = 'plus'
        elif strand=="-":
            color = '255,0,0'
            sName = 'minus'

        return 'track type=wiggle_0 name=%s_%s description=%s_%s color=%s' % (shortname,sName,shortname,sName,color)

    print(wigHeader(shortname,"+"), file=plusHand)
    print(wigHeader(shortname, "-"), file=minusHand)

    for line in handle:
        ref,pos,base,count,reads,quals = line.rstrip().split()
        if ref!=preRef:
            preRef = ref
            print("variableStep chrom=%s" % (ref), file=plusHand)
            print("variableStep chrom=%s" % (ref), file=minusHand)
        if reads.count(".")>0:
            print("%d\t%d" % (int(pos),reads.count(".")), file=plusHand)
        if reads.count(",")>0:
            print("%d\t%d" % (int(pos),reads.count(",")), file=minusHand)

            continue
    plusHand.close()
    minusHand.close()




def getBitValue(n, p):
    """Return the bit at position p of integer n.

    Extracts a single bit at binary position p (zero-indexed from the
    least-significant bit) of the integer n.

    Args:
        n: A non-negative integer to inspect.
        p: Zero-based bit position (0 = least-significant / rightmost bit).

    Returns:
        1 if the bit at position p is set, 0 otherwise.
    """
    return (n >> p) & 1

def strandFlag(flag):
    """Determine the alignment strand from a SAM bitflag value.

    Inspects bit 4 (0x10) of the SAM FLAG field to determine whether a read
    mapped to the reverse strand.

    Args:
        flag: The integer SAM FLAG value (field 2), or a string
            representation of it.

    Returns:
        '+' if bit 4 is 0 (forward strand), '-' if bit 4 is 1 (reverse
        strand), or '*' for any other value.
    """
    flag = int(flag)
    if getBitValue(flag,4)==0:
        return "+"
    elif getBitValue(flag,4)==1:
        return "-"
    else:
        return "*"
