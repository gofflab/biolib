'''
Created on Oct 25, 2009
Misc tools to get information from a SAM/BAM file...
@author: lgoff
'''
from Alignment import Alignment
import intervallib
import os
import pysam
import array
import numpy
import collections
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from inOut.wiggle import WiggleFileWriter

class SAMAlignment(Alignment):
    """Basic object for SAMstring (extends Alignment class)"""
    def __init__(self,readname,chr,start,end,strand,score,readcount,readsequence,cigar,qualstring):
        Alignment.__init__(self,readname,chr,start,end,strand,score=readcount,readcount = readcount,readsequence=readsequence)
        self.qual = qualstring
        self.cigar = cigar

def SAMReader(fname):
    """Iterator for SAMAlignment records (depricated, use pysam)"""
    handle = open(fname,'r')
    for line in handle:
        aln = parseSAMString(line)
        yield aln.toInterval()   

def parseSAMString(samstring):
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
    """Don't use this...it's lazy and it doesn't feel right"""
    handle = open(fname,'r')
    preRef = ''
    prePos = -1
    prePlus = 0
    preMinus = 0
    
    plusHand = open(outDir+shortname+"_plus.wig",'w')
    minusHand = open(outDir+shortname+"_minus.wig",'w')
    
    def wigHeader(shortname,strand):
        if strand=="+":
            color = '0,0,255'
            sName = 'plus'
        elif strand=="-":
            color = '255,0,0'
            sName = 'minus'
        
        return 'track type=wiggle_0 name=%s_%s description=%s_%s color=%s' % (shortname,sName,shortname,sName,color)
    
    print >>plusHand, wigHeader(shortname,"+")
    print >>minusHand, wigHeader(shortname, "-")
    
    for line in handle:
        ref,pos,base,count,reads,quals = line.rstrip().split()
        if ref!=preRef:
            preRef = ref
            print >>plusHand,"variableStep chrom=%s" % (ref)
            print >>minusHand, "variableStep chrom=%s" % (ref)
        if reads.count(".")>0:
            print >>plusHand, "%d\t%d" % (int(pos),reads.count("."))
        if reads.count(",")>0:
            print >>minusHand, "%d\t%d" % (int(pos),reads.count(","))
        
            continue
    plusHand.close()
    minusHand.close()

class Counter:
    """Use in callback function to store read counts within an alignment (includes those that
    are not completely contained within the alignment"""
    mCounts = 0
    def __call__(self,alignment):
        self.mCounts += 1
        
class StrandCounter:
    """Provides a strand-specific number of reads as opposed to total read density"""
    plusCount = 0
    minusCount = 0
    def __call__(self,alignment):
        if alignment.is_reverse:
            self.minusCount += 1
        else:
            self.plusCount += 1


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

def makeCigar():
    pass

def samScanByStrand(samFetch,strand):
    """Generator to iterate over a samFetch using only one of the strands.
    strand should be one of ["+","-"]
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
    if samRead.is_reverse:
        strand = "-"
    else:
        strand = "+"
    return intervallib.Interval(samRead.rname,samRead.pos+1,samRead.pos+1+len(samRead.seq),strand,readcount=1)


def samReadsIntersect(a,b,useStrand = True,offset=0):
    """Checks to see if two samReads (a,b) intersect"""
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
    current = samFetch.next()
    currentInterval = sam2Interval(current)
    
    for x in samFetch:
        next = samFetch.next()
        if samReadsIntersect(current,next,useStrand,offset):
            currentInterval.end = max(currentInterval.end,next.pos+len(next.seq)+1)
            currentInterval.readcount += 1
        else:
            yield currentInterval
            current = samFetch.next()
            currentInterval = sam2Interval(current)    
"""            
def makeContiguousIntervalsByStrand(samHandle,offset=0):
    for strand in ["+","-"]:
        samFetch = samScanByStrand(samHandle.fetch(),strand)
        current = samFetch.next()
        currentInterval = sam2Interval(current)
        
        for next in samFetch:
            if samReadsIntersect(current,next,offset=offset):
                currentInterval.end = max(currentInterval.end,next.pos+len(next.seq)+1)
                currentInterval.readcount += 1
            else:
                yield currentInterval
                current = samFetch.next()
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
    '''
    don't use this function with RNA-seq data because it does not pileup spliced reads properly
    '''
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
    #logger = logging.getLogger(__name__)    
    bamfile = pysam.Samfile(inbamfile, 'rb')

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
        for read_length, count in read_lengths.iteritems():
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
    """This does not work with gapped alignments"""
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
    """This does not work with gapped alignments"""
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
    """This does not work with gapped alignments"""
    
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
            for i in xrange(run):
                cigarMask.append(incrementTable[type])
    return cigarMask

def makePysamCigarMask(cigarTuple,increment=1):
    lookupTable = ['M','I','D','N','S','H','P']
    incrementTable = {
                      'M':increment,
                      'N':0
                      }
    incrementTypes = incrementTable.keys()
    cigarMask = []
    for operation,run in cigarTuple:
        if lookupTable[operation] in incrementTypes:
            for i in xrange(run):
                cigarMask.append(incrementTable[lookupTable[operation]])
    return cigarMask

def bamFetchGappedInterval(bamHandle,chr,start,end,span=1):
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
    """Fetches read density across an interval and finds the start and end position (start and end offset by an index)
     of the kmer with the largest value. Has not been tested yet"""
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
    for i in xrange(end-start+1-k):
        slice = myArr[i:i+k]
        if sum(slice)>maxVal:
            maxVal = sum(slice)
            maxPos = i
    return start+maxPos,end+maxPos

def plotInterval(bamFiles,chr,start,end,name="",span=1,pdfName = "",sumStrands=False):
    nplots = len(bamFiles)
    
    #Setup plot environment
    if not pdfName == "":
        print "Printing figure to %s..." % (pdfName)
        robjects.r.pdf(pdfName,width=8,height=12)
    robjects.r.par(mfrow=array.array('i',[nplots,1]),mar=array.array('i',[2,2,1,0]))
    xaxt = "n"
    count = 0
    for bamFile in bamFiles:
        count+=1
        if count == nplots:
            xaxt = "s"
        baseFname = bamFile.rstrip(".bam")
        bamHandle = pysam.Samfile(bamFile,'rb')
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
    rtrn ={}
    #Fetch total reads in Bam by chromosome
    samfile = pysam.Samfile(bamFile,'rb')
    iter = samfile.fetch(until_eof=True)
    rtrn['readDist'] = {}
    for i in iter:
        rtrn['readDist'][i.rname] = 1 + rtrn['readDist'].get(i.rname,0)
    return rtrn

def getrRNAReads(bamFile,rRNABedFile):
    """Takes a bed file of rRNA genes and queries the bam file to determine the number of unique reads that are mapping to rRNA genes in a given sample"""
    reads = []
    bedIter = intervallib.parseBed(rRNABedFile)
    samfile = pysam.Samfile(bamFile,'rb')
    for bed in bedIter:
        #print "%s\t%s:%d-%d" % (bed.name,bed.chr,bed.start,bed.end)
        res = samfile.fetch(bed.chr,bed.start,bed.end)
        for read in res:
            reads.append(read.qname)
    print "Collapsing to unique"
    return len(uniqify(reads))

def uniqify(seq): 
    # Not order preserving 
    keys = {} 
    for e in seq: 
        keys[e] = 1 
    return keys.keys()

def collapseMatrix(fname):
    """Specifically finds a vector of sums for a chromatin matrix by position"""
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
        vals = numpy.array(map(float,vals))
        sums += vals
        print name
    return names,sums