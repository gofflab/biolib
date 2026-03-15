'''
Miscellaneous tools to get information from a SAM/BAM file.

Provides utilities for parsing SAM/BAM alignment files, computing read
pileups, fetching strand-specific coverage arrays, and plotting read
density across genomic intervals. Built on top of pysam.

Created on Oct 25, 2009

@author: lgoff
'''
import array
import collections
import os

import numpy
import pysam
import rpy2.robjects as robjects

from . import intervallib
from .Alignment import Alignment

# from inOut.wiggle import WiggleFileWriter  # NOTE: inOut.wiggle module not available; WiggleFileWriter commented out

class SAMAlignment(Alignment):
    """Basic object representing a single SAM alignment record.

    Extends the Alignment base class with SAM-specific fields for the
    CIGAR string and base-quality string.

    Attributes:
        qual: Base-quality string from SAM field 11.
        cigar: CIGAR string from SAM field 6 describing the alignment.
    """

    def __init__(self,readname,chr,start,end,strand,score,readcount,readsequence,cigar,qualstring):
        """Initialises a SAMAlignment.

        Args:
            readname: Query template name (SAM field 1).
            chr: Reference sequence name / chromosome (SAM field 3).
            start: 1-based leftmost mapping position (SAM field 4).
            end: Computed end position (start + read length - 1).
            strand: Strand of the alignment, one of '+' or '-'.
            score: Mapping quality score (SAM field 5).
            readcount: Number of reads represented by this alignment
                (typically 1 for a single record).
            readsequence: Read sequence bases (SAM field 10).
            cigar: CIGAR string describing alignment operations (SAM field 6).
            qualstring: ASCII-encoded base-quality string (SAM field 11).
        """
        Alignment.__init__(self,readname,chr,start,end,strand,score=readcount,readcount = readcount,readsequence=readsequence)
        self.qual = qualstring
        self.cigar = cigar

def SAMReader(fname):
    """Iterate over SAM alignment records from a file.

    Deprecated — use pysam directly for new code.

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

    Reads are assumed to be non-paired and non-spliced; the end position is
    derived from the start position plus the read-sequence length.

    Args:
        samstring: A single tab-delimited SAM record line (no trailing
            newline required — it is stripped internally).

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

def pileup2wig(fname,shortname,outDir=os.getcwd()+"/"):
    """Convert a samtools pileup file to strand-specific wiggle files.

    Reads a samtools pileup output file and writes two variableStep wiggle
    files: one for the plus strand (forward reads, indicated by '.') and one
    for the minus strand (reverse reads, indicated by ','). This
    implementation is noted as incomplete / not recommended for production
    use.

    Args:
        fname: Path to the samtools pileup file to read.
        shortname: Base name used for both the wiggle track labels and the
            output file names (``<shortname>_plus.wig`` and
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
            shortname: Base name used in the track name and description fields.
            strand: Strand of the track, either '+' (blue) or '-' (red).

        Returns:
            A wiggle track header string suitable for writing as the first
            line of a wiggle file.
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

class Counter:
    """Callable that accumulates a total read count for use as a pysam callback.

    Designed to be passed as a callback to pysam fetch/pileup methods.
    Counts all reads that overlap the queried region, including those not
    completely contained within it.

    Attributes:
        mCounts: Running total of reads seen so far.
    """

    mCounts = 0

    def __call__(self,alignment):
        """Increment the read counter by one for each alignment seen.

        Args:
            alignment: A pysam AlignedSegment (or compatible) object.
                The alignment itself is not inspected; only its presence
                increments the counter.
        """
        self.mCounts += 1

class StrandCounter:
    """Callable that accumulates strand-specific read counts for use as a pysam callback.

    Separates reads into forward (plus) and reverse (minus) strand tallies
    rather than combining them into a single total.

    Attributes:
        plusCount: Running total of forward-strand reads seen.
        minusCount: Running total of reverse-strand reads seen.
    """

    plusCount = 0
    minusCount = 0

    def __call__(self,alignment):
        """Increment the appropriate strand counter for each alignment seen.

        Args:
            alignment: A pysam AlignedSegment (or compatible) object.
                Strand is determined from the ``is_reverse`` flag.
        """
        if alignment.is_reverse:
            self.minusCount += 1
        else:
            self.plusCount += 1


def getBitValue(n, p):
    """Return the bit at position p of integer n.

    Extracts the single bit at binary position p (zero-indexed from the
    least-significant bit) of the denary integer n.

    Args:
        n: A non-negative integer whose bit is to be inspected.
        p: Zero-based bit position (0 = least-significant / rightmost bit).

    Returns:
        1 if the bit at position p is set, 0 otherwise.
    """
    return (n >> p) & 1

def strandFlag(flag):
    """Determine the alignment strand from a SAM bitflag value.

    Inspects bit 4 (0x10) of the SAM FLAG field to determine whether the
    read mapped to the reverse strand.

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

def makeCigar():
    """Placeholder for CIGAR string construction.

    Not yet implemented.
    """
    pass

def samScanByStrand(samFetch,strand):
    """Yield only reads that map to the specified strand from a pysam fetch iterator.

    Args:
        samFetch: An iterable of pysam AlignedSegment objects, typically
            returned by ``pysam.AlignmentFile.fetch()``.
        strand: The strand to retain.  Must be '+' (forward, non-reverse
            reads) or '-' (reverse reads).

    Yields:
        pysam AlignedSegment objects whose strand matches the requested
        strand value.
    """
    for read in samFetch:
        if strand == "+":
            if read.is_reverse:
                continue
            else:
                yield read
        else:
            if read.is_reverse:
                yield read
            else:
                continue

def sam2Interval(samRead):
    """Convert a pysam AlignedSegment to an intervallib Interval object.

    The interval uses 1-based coordinates (pysam's 0-based ``pos`` is
    incremented by 1) and a readcount of 1.

    Args:
        samRead: A pysam AlignedSegment object with valid ``rname``,
            ``pos``, ``seq``, and ``is_reverse`` attributes.

    Returns:
        An intervallib.Interval representing the read's mapped region,
        with strand set to '+' or '-' according to ``samRead.is_reverse``.
    """
    if samRead.is_reverse:
        strand = "-"
    else:
        strand = "+"
    return intervallib.Interval(samRead.rname,samRead.pos+1,samRead.pos+1+len(samRead.seq),strand,readcount=1)


def samReadsIntersect(a,b,useStrand = True,offset=0):
    """Determine whether two pysam AlignedSegment reads overlap each other.

    Two reads are considered to intersect if their mapped positions overlap
    (allowing for an optional extension by ``offset`` bases). When
    ``useStrand`` is True, reads on different strands or different reference
    sequences are never considered to intersect.

    Args:
        a: A pysam AlignedSegment object.
        b: A pysam AlignedSegment object to compare against ``a``.
        useStrand: If True (default), reads must be on the same reference
            sequence and the same strand (``is_reverse`` must match) to
            be considered intersecting.
        offset: Number of extra bases by which each read's length is
            extended before testing for overlap.  Defaults to 0.

    Returns:
        True if reads a and b overlap (subject to strand and offset rules),
        False otherwise.
    """
    if useStrand:
        if a.rname == b.rname and a.is_reverse == b.is_reverse:
            return not(a.pos>b.pos+len(b.seq)+offset or b.pos>a.pos+len(a.seq)+offset)
        else:
            return False
    else:
        if a.rname == b.rname:
            return not(a.pos>b.pos+len(b.seq)+offset or b.pos>a.pos+len(a.seq)+offset)
        else:
            return False

"""
def makeContiguousIntervals2(samHandle,start='start',end='end',offset=0,useStrand=False):
    '''Generator function to build and iterate over contiguous intervals from a sorted SAM/BAM file.
    If useStrand is True then the function will iterate over one strand at a time.
    '''
    samFetch = samHandle.fetch()
    current = next(samFetch)
    currentInterval = sam2Interval(current)

    for x in samFetch:
        next = next(samFetch)
        if samReadsIntersect(current,next,useStrand,offset):
            currentInterval.end = max(currentInterval.end,next.pos+len(next.seq)+1)
            currentInterval.readcount += 1
        else:
            yield currentInterval
            current = next(samFetch)
            currentInterval = sam2Interval(current)
"""
def makeContiguousIntervalsByStrand(samHandle,offset=0):
    """Generate contiguous genomic intervals from a sorted BAM file, separately per strand.

    Iterates over all reads in the BAM file and merges overlapping reads
    (with optional extension by ``offset``) into contiguous intervals.
    Processing is performed independently for the forward ('+') and reverse
    ('-') strands.

    Args:
        samHandle: An open pysam AlignmentFile object (must be sorted).
        offset: Number of bases by which read extents are extended when
            testing for overlap.  Defaults to 0.

    Yields:
        intervallib.Interval objects representing contiguous merged regions,
        with ``readcount`` reflecting the number of constituent reads.
    """
    for strand in ["+","-"]:
        samFetch = samScanByStrand(samHandle.fetch(),strand)
        current = next(samFetch)
        currentInterval = sam2Interval(current)

        for nxt in samFetch:
            if samReadsIntersect(current, nxt, offset=offset):
                currentInterval.end = max(currentInterval.end, nxt.pos + len(nxt.seq) + 1)
                currentInterval.readcount += 1
            else:
                yield currentInterval
                current = next(samFetch)
                currentInterval = sam2Interval(current)
        yield currentInterval

def generate_pileup_chunks(read_iterator,
                           start, end,
                           unique_only=True,
                           merge_strands=False,
                           fragment_length=-1,
                           dtype=numpy.uint32,
                           max_rlen=2048,
                           chunk_size=8192):
    """Generate read-pileup data in contiguous chunks across a genomic region.

    Iterates over a sorted stream of reads and accumulates per-base read
    depth in fixed-size chunks, yielding each chunk as it is complete.
    Reverse-strand reads may optionally be shifted upstream so that their
    5' end corresponds to the inferred fragment start.

    Note: Do not use with RNA-seq data — spliced reads are not handled
    correctly.

    Args:
        read_iterator: An iterable of pysam AlignedSegment objects sorted
            by position.
        start: 0-based start of the region to pileup.
        end: 0-based (exclusive) end of the region to pileup.
        unique_only: If True (default), reads flagged as PCR/optical
            duplicates (``is_duplicate``) are skipped.
        merge_strands: If True, reverse-strand reads are shifted left by
            ``(read_length - fragment_length)`` bases so both strands
            contribute to the same inferred fragment positions.
        fragment_length: Expected DNA fragment length used to extend reads.
            A value <= 0 means use the actual read length unchanged.
        dtype: numpy dtype for the internal accumulation array.
            Defaults to numpy.uint32.
        max_rlen: Maximum anticipated read length in bases.  The internal
            buffer is sized to accommodate this.  Defaults to 2048.
        chunk_size: Number of bases covered by each yielded chunk.
            Must be >= max_rlen.  Defaults to 8192.

    Yields:
        Tuples of (chunk_start, chunk_end, chunk_array) where chunk_start
        and chunk_end are offsets relative to ``start``, and chunk_array is
        a numpy array of length (chunk_end - chunk_start) containing the
        per-base read depth.
    """
    assert chunk_size >= max_rlen
    assert end > start
    # figure out the boundaries of the first chunk
    chunk_bounds = (start,
                    min(start + chunk_size, end))
    # allocate an array to store the largest possible chunk
    chunk_data = numpy.zeros((chunk_bounds[1] - chunk_bounds[0] + max_rlen,), dtype=dtype)
    chunk_dirty = False
    # reads must be sorted, so implement a check for this
    prev_read_pos = -1
    # iterate through reads
    for read in read_iterator:
        # ignore duplicate reads
        if unique_only and read.is_duplicate:
            continue
        # get attributes from AlignedRead object
        read_start = read.pos
        read_length = read.rlen
        # ensure reads are sorted
        assert read_start >= prev_read_pos
        prev_read_pos = read_start
        # ensure reads aren't too long
        assert read_length <= max_rlen
        # adjust fragment length if too small (can't be smaller than read length)
        if fragment_length <= 0:
            fragment_length = read_length
        # shift the reverse strand reads if the merge_strands option is enabled
        if merge_strands is True:
            if read.is_reverse:
                read_start = max(0, read_start + read_length - fragment_length)
        # now that negative strand tags are shifted, modify the effective read
        # length to the user specified a DNA fragment length
        read_length = fragment_length
        # only consider reads that align within the desired region
        if read_start >= end:
            break
        if (read_start + read_length) > start:
            # if the read starts after the end of the current chunk, need to write the
            # chunk and shift to the next chunk
            while read_start >= chunk_bounds[1]:
                if chunk_dirty:
                    # yield chunk
                    yield (chunk_bounds[0]-start), (chunk_bounds[1]-start), chunk_data[0:(chunk_bounds[1]-chunk_bounds[0])]
                    # add chunk to array
                    #arr[chunk_bounds[0]-start:chunk_bounds[1]-start] += chunk_data[0:(chunk_bounds[1]-chunk_bounds[0])]
                    # shift end of chunk to beginning of next chunk and clear rest of array
                    chunk_data[0:max_rlen] = chunk_data[-max_rlen:]
                    chunk_data[max_rlen:] = 0
                    # check if chunk no longer dirty
                    chunk_dirty = chunk_data[0:max_rlen].any()
                # get next chunk
                chunk_bounds = (chunk_bounds[0] + chunk_size,
                                min(chunk_bounds[1] + chunk_size, end))
            # add coverage from the current read
            chunk_data[max(0, read_start - chunk_bounds[0]):read_start + read_length - chunk_bounds[0]] += 1
            chunk_dirty = True
    # flush last chunks
    while chunk_dirty and (chunk_bounds[0] < end):
        # yield chunk
        yield (chunk_bounds[0]-start), (chunk_bounds[1]-start), chunk_data[0:(chunk_bounds[1]-chunk_bounds[0])]
        # shift end of chunk to beginning of next chunk and clear rest of array
        chunk_data[0:max_rlen] = chunk_data[-max_rlen:]
        chunk_data[max_rlen:] = 0
        # check if chunk no longer dirty
        chunk_dirty = chunk_data[0:max_rlen].any()
        # get next chunk
        chunk_bounds = (chunk_bounds[0] + chunk_size,
                        min(chunk_bounds[1] + chunk_size, end))
    # delete chunk array
    del chunk_data


def bam_to_wiggle(inbamfile, wigfile,
                  unique_only=False,
                  merge_strands=False,
                  fragment_length=-1,
                  norm=False):
    """Convert a BAM file to a compressed wiggle file.

    Computes per-base read depth across every reference sequence in the BAM
    file and writes the result as a wiggle file using WiggleFileWriter (from
    the inOut.wiggle module).  Note: WiggleFileWriter is currently
    unavailable — calling this function will raise a NameError.

    Args:
        inbamfile: Path to the input BAM file (must be sorted and indexed).
        wigfile: Path to the output wiggle file to write.
        unique_only: If True, reads flagged as PCR/optical duplicates are
            excluded from the pileup.  Defaults to False.
        merge_strands: If True, reverse-strand reads are shifted upstream
            so both strands reflect inferred fragment start positions.
            Defaults to False.
        fragment_length: Expected DNA fragment length used to extend reads.
            A value <= 0 means use the actual read length unchanged.
        norm: If True, read depths are normalised to reads-per-kilobase per
            million mapped reads (RPKM-style).  Defaults to False.
    """
    #logger = logging.getLogger(__name__)
    bamfile = pysam.AlignmentFile(inbamfile, 'rb')

    # count reads and get other info from BAM file
    reads = 0
    read_lengths = collections.defaultdict(lambda: 0)
    for read in bamfile.fetch():
        # only consider mapped reads for now
        if read.is_unmapped:
            continue
        reads += 1
        read_lengths[read.rlen] += 1
    # find normalization factor
    if norm == True:
        # find best read length
        best_read_length, best_count = 0, 0
        for read_length, count in read_lengths.items():
            if count > best_count:
                best_count = count
                best_read_length = read_length
        assert best_read_length > 0
        # normalize by read length and number of mapped reads
        norm_factor = (1.0e9) / (max(best_read_length, fragment_length) * reads)
    else:
        norm_factor = 1.0

    refs = bamfile.references
    lengths = bamfile.lengths
    # NOTE: WiggleFileWriter is unavailable (inOut.wiggle not importable); this will raise NameError if called
    wigglewriter = WiggleFileWriter(wigfile, compress=True, span=10)
    # convert each chromosome to wiggle
    for ref, length in zip(refs, lengths):
        # pileup the reads chunks at a time
        for pileupchunk in generate_pileup_chunks(bamfile.fetch(ref),
                                                  start=0,
                                                  # TODO: some wiggle writing error with length going past limit
                                                  end=length - max(0, fragment_length),
                                                  unique_only=unique_only,
                                                  merge_strands=merge_strands,
                                                  fragment_length=fragment_length,
                                                  chunk_size=1048576):
            chunk_start, chunk_end, chunk_data = pileupchunk
            if norm == True:
                chunk_data *= norm_factor
            #wigglewriter.write_variable_step(ref, chunk_start, chunk_end, chunk_data)
            wigglewriter.write_span(ref, chunk_start, chunk_end, chunk_data)
        #logger.debug("BAM %s -> WIG %s chromsome %s finished" % (inbamfile, wigfile, ref))
    # wiggle file done
    wigglewriter.close()
    # done with BAM file
    bamfile.close()

def bamFetchFlank(bamHandle,chr,pos,flankSize=1000,fragment_length=200):
    """Compute merged-strand read-depth in a window centred on a genomic position.

    Fetches reads from a BAM file within ``pos ± (flankSize + fragment_length)``
    and accumulates per-base coverage into a numpy array.  Reverse-strand
    reads are shifted upstream to align with their inferred fragment start.

    Note: Does not handle gapped (spliced) alignments correctly.

    Args:
        bamHandle: An open pysam AlignmentFile object.
        chr: Reference sequence name / chromosome to query.
        pos: Centre position (0-based) of the window.
        flankSize: Number of bases to include on each side of ``pos`` in the
            returned array.  Defaults to 1000.
        fragment_length: Expected DNA fragment length used to extend reverse-
            strand reads.  A value <= 0 means use the actual read length.
            Defaults to 200.

    Returns:
        A numpy array of length ``2 * flankSize + 1`` containing the
        per-base read depth centred on ``pos``.
    """
    #Create container to hold pos +- (flankSize+fragment_length)
    arr = numpy.zeros(2*(flankSize+fragment_length)+1)
    range = (pos-flankSize-fragment_length,pos+flankSize+fragment_length)

    readIter = bamHandle.fetch(chr,range[0],range[1])
    for read in readIter:
        if read.is_unmapped:
            continue
        read_start = read.pos
        read_length = read.rlen

        if fragment_length <= 0:
            fragment_length = read_length
        if read.is_reverse:
            read_start = max(0, read_start + read_length - fragment_length)
        # now that negative strand tags are shifted, modify the effective read
        # length to the user specified a DNA fragment length
        read_length = fragment_length
        # only consider reads that align within the desired region
        arr[max(0, read_start - range[0]):read_start + read_length - range[0]] += 1
    return arr[fragment_length:fragment_length+2*flankSize+1]

def bamFetchFlank_byStrand(bamHandle,chr,pos,flankSize=1000,fragment_length=200,span=1):
    """Compute strand-specific read-depth arrays in a window centred on a genomic position.

    Similar to ``bamFetchFlank`` but returns separate arrays for the sense
    (forward) and antisense (reverse) strands.  Reverse-strand reads are
    extended to the inferred fragment start when ``fragment_length`` exceeds
    the read length.

    Note: Does not handle gapped (spliced) alignments correctly.

    Args:
        bamHandle: An open pysam AlignmentFile object.
        chr: Reference sequence name / chromosome to query.
        pos: Centre position (0-based) of the window.
        flankSize: Number of bases to include on each side of ``pos`` in
            each returned array.  Defaults to 1000.
        fragment_length: Expected DNA fragment length used to extend reverse-
            strand reads.  A value <= 0 means use the actual read length.
            Defaults to 200.
        span: Step size for down-sampling the output arrays.  A value of 1
            (default) returns every base; 2 returns every other base, etc.

    Returns:
        A tuple (senseArr, antisenseArr) where each element is a numpy
        array of length ``(2 * flankSize + 1) / span`` containing per-base
        read depth for the respective strand, centred on ``pos``.
    """
    senseArr = numpy.zeros(2*(flankSize+fragment_length)+1)
    antisenseArr = numpy.zeros(2*(flankSize+fragment_length)+1)

    range = (pos-flankSize-fragment_length,pos+flankSize+fragment_length)


    readIter = bamHandle.fetch(chr,range[0],range[1])
    for read in readIter:
        if read.is_unmapped:
            continue
        read_start = read.pos
        read_length = read.rlen

        if not read.is_reverse:
            if fragment_length <= 0:
                fragment_length = read_length

            read_length = fragment_length
            senseArr[max(0,read_start - range[0]):read_start + read_length - range[0]] += 1
        else:
            read_end = read_start + read_length
            if fragment_length > read_length:
                read_start = max(0,read_start + read_length - fragment_length)
            antisenseArr[max(0,read_start-range[0]):read_end - range[0]] += 1
    return (senseArr[fragment_length:fragment_length+2*flankSize+1:span],antisenseArr[fragment_length:fragment_length+2*flankSize+1:span])

def bamFetchInterval(bamHandle,chr,start,end,fragment_length=200,span=1):
    """Compute strand-specific read-depth arrays across a genomic interval.

    Fetches reads from the BAM file that overlap ``[start, end]`` and
    accumulates per-base read depth separately for the sense and antisense
    strands.  Reverse-strand reads whose actual length is less than
    ``fragment_length`` are extended upstream to the inferred fragment start.

    Note: Does not handle gapped (spliced) alignments correctly.

    Args:
        bamHandle: An open pysam AlignmentFile object.
        chr: Reference sequence name / chromosome to query.
        start: 0-based start of the interval.
        end: 0-based end of the interval (inclusive).
        fragment_length: Expected DNA fragment length used to extend reads.
            A value <= 0 means use the actual read length unchanged.
            Defaults to 200.
        span: Step size for down-sampling the output arrays.  Defaults to 1.

    Returns:
        A tuple (senseArr, antisenseArr) where each element is a numpy
        array of length ``(end - start + 1) / span`` containing per-base
        read depth for the respective strand across the interval.
    """

    senseArr = numpy.zeros(end-start+(2*fragment_length)+1)
    antisenseArr = numpy.zeros(end-start+(2*fragment_length)+1)

    range = (start-fragment_length,end+fragment_length)
    intervalSize = end-start+1

    readIter = bamHandle.fetch(chr,range[0],range[1])
    for read in readIter:
        if read.is_unmapped:
            continue
        read_start = read.pos
        read_length = read.rlen

        if not read.is_reverse:
            if fragment_length <=0:
                fragment_length = read_length

            read_length = fragment_length
            senseArr[max(0,read_start - range[0]):read_start + read_length - range[0]] += 1
        else:
            read_end = read_start+read_length
            if fragment_length > read_length:
                read_start = max(0,read_start + read_length - fragment_length)
            antisenseArr[max(0,read_start-range[0]):read_end - range[0]] += 1
    return(senseArr[fragment_length:fragment_length+intervalSize:span],antisenseArr[fragment_length:fragment_length+intervalSize:span])

def makeCigarMask(cigar,increment=1):
    """Build a per-base mask vector from a CIGAR string.

    Parses a text CIGAR string and produces a flat list where each element
    corresponds to one reference base consumed by the alignment.  'M'
    (match/mismatch) operations contribute ``increment`` to each position;
    'N' (intron/skip) operations contribute 0.  Other CIGAR operations that
    do not consume reference bases (e.g. 'I', 'S', 'H', 'P') are omitted
    from the output.

    Args:
        cigar: A CIGAR string such as ``'36M'`` or ``'20M1000N16M'``.
        increment: Value assigned to each matched ('M') reference base in
            the output mask.  Defaults to 1.

    Returns:
        A list of numeric values (each 0 or ``increment``) with one entry
        per reference base consumed by the alignment.
    """
    incrementTable = {
                      'M':increment,
                      'N':0
                      }
    incrementTypes = incrementTable.keys()
    #Split CIGAR into components
    components = []
    values = ['M','I','D','N','S','H','P']
    tmpString = ''
    for i in cigar:
        #print i
        tmpString += i
        if i in values:
            components.append(tmpString)
            tmpString = ''
    components = [(x[-1],int(x[:-1])) for x in components]
    #make vector to hold full length and increment appropriate indices
    cigarMask = []
    for type,run in components:
        if type in incrementTypes:
            for i in range(run):
                cigarMask.append(incrementTable[type])
    return cigarMask

def makePysamCigarMask(cigarTuple,increment=1):
    """Build a per-base mask vector from a pysam CIGAR tuple.

    Equivalent to ``makeCigarMask`` but accepts the pysam representation
    of a CIGAR string (a list of (operation_code, length) integer pairs)
    rather than a text CIGAR string.  'M' operations contribute
    ``increment``; 'N' operations contribute 0; other operations that do
    not consume reference bases are omitted.

    Args:
        cigarTuple: A sequence of (operation, length) pairs as returned by
            pysam's ``AlignedSegment.cigar`` attribute.  Operation codes
            follow the SAM spec order: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P.
        increment: Value assigned to each matched ('M') reference base.
            Defaults to 1.

    Returns:
        A list of numeric values (each 0 or ``increment``) with one entry
        per reference base consumed by the alignment.
    """
    lookupTable = ['M','I','D','N','S','H','P']
    incrementTable = {
                      'M':increment,
                      'N':0
                      }
    incrementTypes = incrementTable.keys()
    cigarMask = []
    for operation,run in cigarTuple:
        if lookupTable[operation] in incrementTypes:
            for i in range(run):
                cigarMask.append(incrementTable[lookupTable[operation]])
    return cigarMask

def bamFetchGappedInterval(bamHandle,chr,start,end,span=1):
    """Compute strand-specific read-depth arrays across an interval, respecting CIGAR gaps.

    Unlike ``bamFetchInterval``, this function uses each read's CIGAR
    information (via ``makePysamCigarMask``) so that intronic regions ('N'
    operations) do not contribute to the depth.  Fragment-length extension
    is not yet implemented (TODO).

    Args:
        bamHandle: An open pysam AlignmentFile object.
        chr: Reference sequence name / chromosome to query.
        start: 0-based start of the interval.
        end: 0-based end of the interval (inclusive).
        span: Step size for down-sampling the output arrays.  Defaults to 1.

    Returns:
        A tuple (senseArr, antisenseArr) where each element is a numpy
        array of length ``(end - start + 1) / span`` containing per-base
        read depth for the respective strand across the interval.
    """
    #TODO incoporate fragment size into reads (see above), default 200nt
    intervalSize = end-start+1
    senseArr = numpy.zeros(intervalSize)
    antisenseArr = numpy.zeros(intervalSize)

    readIter = bamHandle.fetch(chr,start,end)
    for read in readIter:
        if read.is_unmapped:
            continue
        readStart = read.pos
        #print read.cigar
        mask = makePysamCigarMask(read.cigar,increment=1)
        if not read.is_reverse:
            #This is necessary to avoid a left overhang when merging the cigarmask with the final vector
            leftOffset = -min(0,readStart-start)
            """
            if readStart-start < 0:
                leftOffset = -(readStart-start)
            else:
                leftOffset = 0

            Debugging...

            #print read.pos #(this is the problem Samtools takes reads that start before 'start')
            print readStart-start
            print mask
            print leftOffset
            print len(senseArr[readStart-start:readStart-start+len(mask)])
            print len(mask[:len(senseArr)-(readStart-start)])
            """
            senseArr[readStart-start+leftOffset:readStart-start+len(mask)] += mask[leftOffset:len(senseArr)-(readStart-start)]
        else:
            leftOffset = -min(0,readStart-start)
            antisenseArr[readStart-start+leftOffset:readStart-start+len(mask)] += mask[leftOffset:len(antisenseArr)-(readStart-start)]
    return senseArr[::span],antisenseArr[::span]

def findLargestKmer(bamHandle,chr,start,end,strand,k=21,gapped=False,span=1):
    """Find the k-mer window with the highest total read depth within an interval.

    Computes per-base read depth across the interval (using either the
    simple or gapped pileup function) and slides a window of size ``k``
    across the appropriate strand array to locate the window whose summed
    depth is largest.

    Note: This function has not been tested yet.

    Args:
        bamHandle: An open pysam AlignmentFile object.
        chr: Reference sequence name / chromosome to query.
        start: 0-based start of the interval.
        end: 0-based end of the interval (inclusive).
        strand: Which strand array to search; '+' uses the sense array,
            '-' uses the antisense array.
        k: Window size in bases.  Defaults to 21.
        gapped: If True, uses ``bamFetchGappedInterval`` (CIGAR-aware
            pileup); otherwise uses ``bamFetchInterval``.  Defaults to False.
        span: Down-sampling step passed to the pileup function.
            Defaults to 1.

    Returns:
        A tuple (window_start, window_end) giving the genomic coordinates
        of the highest-scoring k-mer window.  Both values are offset from
        ``start`` by the index of the best window.
    """
    if not gapped:
        sense,antisense = bamFetchInterval(bamHandle,chr,start,end,span=span)
    else:
        sense,antisense = bamFetchGappedInterval(bamHandle,chr,start,end,span=span)

    if strand == "+":
        myArr = sense
    elif strand == "-":
        myArr = antisense

    maxVal = 0
    maxPos = -1
    for i in range(end-start+1-k):
        slice = myArr[i:i+k]
        if sum(slice)>maxVal:
            maxVal = sum(slice)
            maxPos = i
    return start+maxPos,end+maxPos

def plotInterval(bamFiles,chr,start,end,name="",span=1,pdfName = "",sumStrands=False):
    """Plot read depth across a genomic interval for one or more BAM files.

    Uses rpy2 to create a multi-panel line plot, one panel per BAM file.
    Forward-strand depth is shown in blue (positive y-axis) and reverse-
    strand depth in red (negative y-axis) unless ``sumStrands`` is True,
    in which case a single combined black trace is drawn.  Optionally saves
    the plot to a PDF.

    Args:
        bamFiles: A list of paths to BAM files to plot (one panel each).
        chr: Reference sequence name / chromosome to display.
        start: 0-based start of the display window.
        end: 0-based end of the display window (inclusive).
        name: Optional label appended to each panel title.  Defaults to ''.
        span: Down-sampling step passed to the pileup function.
            Defaults to 1.
        pdfName: If non-empty, the plot is written to this PDF path; otherwise
            an interactive R window is opened.  Defaults to ''.
        sumStrands: If False (default), sense and antisense tracks are
            plotted separately with opposite sign.  If True, strand depths
            are summed into a single positive trace.
    """
    nplots = len(bamFiles)

    #Setup plot environment
    if not pdfName == "":
        print("Printing figure to %s..." % (pdfName))
        robjects.r.pdf(pdfName,width=8,height=12)
    robjects.r.par(mfrow=array.array('i',[nplots,1]),mar=array.array('i',[2,2,1,0]))
    xaxt = "n"
    count = 0
    for bamFile in bamFiles:
        count+=1
        if count == nplots:
            xaxt = "s"
        baseFname = bamFile.rstrip(".bam")
        bamHandle = pysam.AlignmentFile(bamFile,'rb')
        sense,antisense = bamFetchGappedInterval(bamHandle,chr,start,end,span=span)

        if sumStrands == False:
            lim = int(max(numpy.concatenate((sense,antisense))))
            robjects.r.plot(range(start,end+1,span),sense,type="l",ylab="",xlab="",main=baseFname+"_"+name,col="blue",ylim=array.array('i',[-lim,lim]),xlim=array.array('i',[start,end]),xaxt=xaxt,bty='n')
            #robjects.r.polygon(range(start,end+1,span),sense,col="blue",border=False)
            robjects.r.lines(range(start,end+1,span),-antisense,col="red")
            robjects.r.abline(h=0,col="black")
        elif sumStrands == True:
            myVec = sense+antisense
            robjects.r.plot(range(start,end+1,span),myVec,type="l",ylab="",xlab="",main=baseFname+"_"+name,col="black",ylim=array.array('i',[0,int(max(myVec))]),xlim=array.array('i',[start,end]),xaxt=xaxt,bty='n')

    if not pdfName == "":
        robjects.r['dev.off']()

def bamStats(bamFile):
    """Compute per-chromosome read counts for a BAM file.

    Iterates over every read in the BAM file (including unmapped reads) and
    tallies how many reads map to each reference sequence.

    Args:
        bamFile: Path to the BAM file.

    Returns:
        A dict with a single key ``'readDist'`` whose value is itself a
        dict mapping reference sequence index (``rname``) to the number of
        reads mapping to that reference.
    """
    rtrn ={}
    #Fetch total reads in Bam by chromosome
    samfile = pysam.AlignmentFile(bamFile,'rb')
    iter = samfile.fetch(until_eof=True)
    rtrn['readDist'] = {}
    for i in iter:
        rtrn['readDist'][i.rname] = 1 + rtrn['readDist'].get(i.rname,0)
    return rtrn

def getrRNAReads(bamFile,rRNABedFile):
    """Count unique reads that map to rRNA gene loci.

    Parses a BED file of rRNA gene coordinates and queries the BAM file for
    each locus, collecting all overlapping read names.  Duplicate read names
    are collapsed before returning the final count.

    Args:
        bamFile: Path to the sorted, indexed BAM file to query.
        rRNABedFile: Path to a BED file listing rRNA gene intervals.

    Returns:
        The number of unique read names (query names) that overlap at least
        one rRNA gene locus.
    """
    reads = []
    bedIter = intervallib.parseBed(rRNABedFile)
    samfile = pysam.AlignmentFile(bamFile,'rb')
    for bed in bedIter:
        #print "%s\t%s:%d-%d" % (bed.name,bed.chr,bed.start,bed.end)
        res = samfile.fetch(bed.chr,bed.start,bed.end)
        for read in res:
            reads.append(read.qname)
    print("Collapsing to unique")
    return len(uniqify(reads))

def uniqify(seq):
    """Return the unique elements of a sequence (order not preserved).

    Args:
        seq: Any iterable of hashable elements.

    Returns:
        A dict_keys view containing one entry per unique element found in
        ``seq``.  The original order is not preserved.
    """
    # Not order preserving
    keys = {}
    for e in seq:
        keys[e] = 1
    return keys.keys()

def collapseMatrix(fname):
    """Sum a tab-delimited chromatin matrix column-wise across all samples.

    Reads a matrix file whose first row is a header and whose subsequent
    rows each begin with two identifier fields (sample and name) followed
    by numeric values.  Returns the element-wise sum of all data rows and
    the list of row names.

    Args:
        fname: Path to a tab-delimited matrix file.  Expected format: the
            first line is a header whose columns (after the leading
            identifier columns) name the positions.  Each subsequent line
            starts with a sample identifier and a row name, followed by
            numeric values.

    Returns:
        A tuple (names, sums) where ``names`` is a list of row-name strings
        (second column of each data row) and ``sums`` is a numpy array of
        the column-wise sums across all data rows.
    """
    handle = open(fname,'r')
    header = handle.readline().rstrip()
    header = header.split("\t")[1:]
    sums = numpy.zeros(len(header))
    names = []

    for line in handle:
        vals = line.rstrip().split("\t")
        sample = vals.pop(0)
        name = vals.pop(0)
        names.append(name)
        vals = numpy.array([float(x) for x in vals])
        sums += vals
        print(name)
    return names,sums
