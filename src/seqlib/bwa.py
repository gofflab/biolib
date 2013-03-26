'''
Created on Jul 30, 2009
Python wrappers for BWA algorithm

@author: lgoff

Example of base commands:
BWA Align:
    bwa aln -c /seq/compbio-hp/lgoff/genomes/hg18/hg18.fa test.fastq >test.sai
BWA SAMSE:
     bwa samse /seq/compbio-hp/lgoff/genomes/hg18/hg18.fa test.sai test.fastq
'''
import os,copy
from Alignment import *

prefix = "/seq/compbio-hp/lgoff/genomes/hg18/hg18.fa"
ref_index = prefix+".fai"

#=================
class SAMAlignment(Alignment):
    def __init__(self,readname,chr,start,end,strand,score,readcount,readsequence,cigar,qualstring):
        Alignment.__init__(self,readname,chr,start,end,strand,score=readcount,readcount = readcount,readsequence=readsequence)
        self.qual = qualstring
        self.cigar = cigar

def SAMReader(fname):
    """Iterator for SAMAlignment records"""
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

def joinSAMIntervals(iter,start='start',end='end',offset=0):
    """
    Returns a list of independent non-overlapping intervals for each strand overlapping by offset if set
    ***SAM file must first be sorted using 'samtools sort'***
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
        print i
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
    for fname in files:
        shortname = fname.rstrip(".fastq")
        command = "bsub -q %s -N -o /dev/null -P BWA_Align 'bwa aln -c -n %d %s %s >%s.sai 2>%s.e'" % (queue,mismatches,prefix,fname,shortname,shortname)
        os.system(command)
    return

def bwaSamseSubmit(files,mismatches=2,queue='broad'):
    for fname in files:
        shortname = fname.rstrip(".sai")
        command = "bsub -q %s -N -o /dev/null -P BWA_Samse 'bwa samse %s %s.sai %s.fastq >%s.sam 2>%s.e'" % (queue,prefix,shortname,shortname,shortname,shortname)
        os.system(command)
    return

def makeBam(files,queue='broad'):
    for fname in files:
        shortname = fname.rstrip("*.sam")
        command = "bsub -q %s -N -o /dev/null -P SAM2BAM 'samtools view -h -bt %s -o %s.bam %s 2>%s.bam.e'" % (queue,ref_index,shortname,fname,shortname)
        os.system(command)
    return

def samSort(files,queue='broad'):
    for fname in files:
        shortname = fname.rstrip("*.bam")+"_sorted"
        command = "samtools sort %s %s" % (fname,shortname)
        print "Sorting file: %s" % fname
        os.system(command)
    return



def pileup2wig(fname,shortname,outDir=os.getcwd()+"/"):
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