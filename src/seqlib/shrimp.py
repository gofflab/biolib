#!/usr/bin/python
import string,os,random,sys,glob,solid
from subprocess import *
from intervallib import *
from Alignment import *
import genomelib

###############
#SHRiMP Program Variables
#
###############
shrimp_dir = '/home/radon01/lgoff/apps/SHRiMP/bin/'
utils_dir = '/home/radon01/lgoff/apps/SHRiMP/utils/'
app = 'rmapper-cs'
mismatches = 2
max_hits = 10
max_length = 25

order = ["readname","contigname","strand","contigstart","contigend","readstart","readend","readlength","score","editstring","readsequence"]
#######################
class ShrimpRead(Alignment):
    """Extends Alignment class to include a few SHRiMP-specific attributes and methods"""
    
    def __init__(self,readname,chr,start,end,strand,readstart,readend,score,readcount,readsequence,editstring,readlength):
        Alignment.__init__(self,readname,chr,start,end,strand,score,readcount,readsequence)
        self.readstart = int(readstart)
        self.readend = int(readend)
        self.readcount = int(readcount)
        self.editstring = editstring
        self.readlength = readlength
        self.crossovers = self.editstring.count("x")
        self.nSNPs = 0
        for s in ['A','C','G','T']:
            self.nSNPs += self.editstring.count(s)
        self.aligner = "shrimp"
        
    def __len__(self):
        return self.readlength
    
    def __str__(self):
        return "SHRiMP:%s:%s:%d:%d" % (self.readname,self.chr,self.start,self.end)
    
    def shrimpString(self):
        return ">%s_x%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" % (self.readsequence,self.readcount,self.chr,self.strand,self.start,self.end,self.readstart,self.readend,self.readlength,self.score,self.editstring,self.readsequence)

class ProbCalcRead(ShrimpRead):
    """Extends ShrimpRead class to include statistical output from probcalc"""
    def __init__(self,readname,chr,start,end,strand,readstart,readend,score,readcount,editstring,readlength,normodds,pgenome,pchance,readsequence=''):
        ShrimpRead.__init__(self,readname,chr,start,end,strand,readstart,readend,score,readcount,readsequence,editstring,readlength)
        self.readsequence = misc.nuID2seq(self.readname)
        self.normodds = float(normodds)
        self.pgenome = float(pgenome)
        self.pchance = float(pchance)
        
def prepShrimp(file,basedir,binSize=1000):
    """Takes as input a .csfasta file (only) and optionally the number of reads per output file (binSize).
    Takes input .csfasta file and creates subdirectories for reads and results.  Then parses the .csfasta
    file into files of 'binSize' reads.
    """
        
    if not file.endswith('.csfasta'):
        raise ValueError("prepShrimp requires a .csfasta file as input")
    #sample = file.rstrip('.csfasta')
    readsDir = basedir+"/reads/"
    os.mkdir(readsDir)
    resultsDir = basedir+"/results/"
    os.mkdir(resultsDir)
    curDir=os.getcwd()
    os.chdir(readsDir)
    command = "python %ssplitreads.py %d %s" % (utils_dir,binSize,curDir+"/"+file)
    os.system(command)
    os.chdir(curDir)

def GenRandom(length = 10, chars=string.letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])

def submitShrimp(queue="broad",cwd = os.getcwd(),outDir="../results/",readLength=25):
    """Takes all .csfasta files in cwd and submits them as jobs to LSF.  Assumes current directory is
    'reads' and results are in '../results'.  Tags stdout as 'file.o' and stderr as 'file.e'
    """
    files = os.listdir(cwd)
    for file in files:
        if file.endswith('.csfasta'):
            basename = file.rstrip('.csfasta')
            call = """bsub -P "SHRiMP" -q %s -o /dev/null "%s -r %d -o %d -R %s %s >%s 2>%s" """ % (queue,shrimp_dir+app,readLength,max_hits,file,genomelib.genome_dir+"/"+genomelib.genome_file,outDir+basename+'.o',outDir+basename+'.e')
            os.system(call)

def parseShrimp(handle):
    """
    This is the second version to parse SHRiMP v1.1 output or later. (Added 9-19-08)
    Generator Function to return iterator of ShrimpRead objects SHRiMP alignment output file
    """
    
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break
    while True:
        if line[0] <>">":
            raise ValueError("Records in Fasta files should start with a '>' character")
        #Split row into list
        parsedList = line[1:].rstrip().split("\t")
        read = ShrimpRead(parsedList[10],
                          parsedList[1],
                          parsedList[3],
                          parsedList[4],
                          parsedList[2],
                          parsedList[5],
                          parsedList[6],
                          parsedList[8],
                          parsedList[0].split("_x")[1],
                          parsedList[10],
                          parsedList[9],
                          parsedList[7]
                          )
        
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            line = handle.readline()
        
        yield read
        if not line : return #StopIteration
    assert False, "Should not reach this line"

def parseProbcalc(handle):
    """
    This is the second version to parse SHRiMP v1.1 PROBCALC output or later. (Added 9-11-09)
    Generator Function to return iterator of ShrimpRead objects from probcalc output file
    """
    
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break
    while True:
        if line[0] <>">":
            raise ValueError("Records in Fasta files should start with a '>' character")
        #Split row into list
        parsedList = line[1:].rstrip().split("\t")
        read = ProbCalcRead(parsedList[0].split("_x")[0],
                          parsedList[1],
                          parsedList[3],
                          parsedList[4],
                          parsedList[2],
                          parsedList[5],
                          parsedList[6],
                          parsedList[8],
                          parsedList[0].split("_x")[1],
                          parsedList[9],
                          parsedList[7],
                          parsedList[10],
                          parsedList[11],
                          parsedList[12]
                          )
        
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            line = handle.readline()
        
        yield read
        if not line : return #StopIteration
    assert False, "Should not reach this line"

'''
def joinShrimp():
    """Deprecated"""
    
    import glob
    files = glob.glob("*.aln")
    
    strands = ['+','-']
    
    for file in files:
        handle = open(file,'r')
        id = file.rstrip(".aln")
        sampleName,chrom = id.split("_")
        outputs = {}
        for strand in strands:
            outputs[strand] = open(sampleName+"_"+chrom+"_"+strand+".bed",'w')
    
        intervals_plus = []
        intervals_minus = []
        
        keys = ['chr','start','end','label','count','strand']
        for line in handle:
            if line.startswith("#"):
                continue
            values = line.rstrip().split("\t")
            result = dict(zip(keys,values))
            (result['start'],result['end'],result['count']) = (int(result['start']),int(result['end']),int(result['count']))
            if result['strand'] == '+':
                intervals_plus.append(result)
            if result['strand'] == '-':
                intervals_minus.append(result)
        
        intervals_plus = genome.join_intervals_sum(intervals_plus)
        intervals_minus = genome.join_intervals_sum(intervals_minus)
        
        for interval_plus in intervals_plus:
            print >> outputs['+'], "%s\t%d\t%d\t.\t%d\t+" % (chrom,interval_plus['start'],interval_plus['end'],interval_plus['count'])
        for interval_minus in intervals_minus:
            print >> outputs['-'], "%s\t%d\t%d\t.\t%d\t-" % (chrom,interval_minus['start'],interval_minus['end'],interval_minus['count'])
'''

def split_shrimp(fname,outDir=os.getcwd()):
    """Takes the concatenated output of SHRiMP outputs (fname) and splits it into multiple files by number of mismatches per alignment """
    handle = open(fname,'r')
    iter = parseShrimp(handle)
    for i in iter:
        h2=open(outDir+"/"+fname+'_'+str(i.crossovers)+".split",'a')
        h2.write(i.shrimpString())
        h2.close
    return

