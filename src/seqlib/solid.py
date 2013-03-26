#!/usr/bin/python
import sys,os
#import math
import misc
#from random import choice
#import string

#######
#Generic SOLiD parameters
########
P2_seq = 'CGCCTTGGCCGTACAGCA'
P2_CS_seq = 'C330201030313112312'

mapping = {}
mapping["0"] = {"T":"T","A":"A","C":"C","G":"G"}
mapping["1"] = {"T":"G","A":"C","C":"A","G":"T"}
mapping["2"] = {"T":"C","A":"G","C":"T","G":"A"}
mapping["3"] = {"T":"A","A":"T","C":"G","G":"C"}

def linker_oligos(linker = P2_seq):
    '''returns a set of oligo-nts that constitute SOLiD P2 adapter sequence'''
    oligos = []
    for i in range(len(linker)):
        oligos.append(linker[:i+1])
    return set(oligos)

#################################################################
#CSSeq Class definition:  Basic class of Colorspace sequence
#################################################################
class CSSeq:
    "Defines the basic sequence class for the pipeline (DNA or CS)"
    def __init__(self,name,sequence,readcount=1):
        self.name = name
        self.sequence = sequence
        self.readcount = readcount
        self.matches = []
        self.qual = []
        self.space = "CS"
        self.trimmed = False
        #self.count = 0 
    
    def __len__(self):
        return len(self.sequence)
    
    def __str__(self):
        return self.sequence
    def __repr__(self):
        return self.name
    
#    def __repr__(self):
#        return "***Object of class 'CSSeq'***\nName:     %s\nSequence: %s\nSpace:    %s\nTrimmed:  %s" % (self.name,self.sequence,self.space,self.trimmed)
    
    #Added per request by Ron to add IVGN samples to database from .csfasta
    #def SQLOutput(self):
    #    """Returns string of BeadName<tab>CSsequence<tab>DNAsequence for insert into database"""
    #    CSseq = self.sequence
    #    self.CSToDNA()
    #    DNAseq = self.sequence
    #    return ('%s\t%s\t%s\t' % (self.name,CSseq,self.sequence))
        
    def returnFasta(self):
        return ('>%s\n%s' % (self.name,self.sequence))
    
    def returnSHRiMPcsfasta(self):
        return ('>%s_x%d\n%s') % (self.name,self.readcount,self.sequence)
    
    def returnQual(self):
        return('>%s\n%s' % (self.name," ".join(q for q in self.qual)))
    
    def printFasta(self):
        print ('>%s\n%s' % (self.name,self.sequence))
    
    def CSToDNA(self):
        """
        This function will convert the colorspace 'self.sequence' to DNA space
        """
        if self.space!="CS":
            raise TypeError('Not a colorspace sequence')
        
        res = ''
        letter = ''
        
        for i in self.sequence:
            if (letter == ''):
                letter = res = i
                continue
            else:
                letter = mapping[i][letter]
            res = res + letter
        self.space = "DNA"
        self.sequence = res[1:]

    def strip_solid_linker(self, linker=None):
        '''strips the linker from a SOLiD read in DNA-space or colorspace'''
        if self.trimmed:
            raise TypeError('Sequence has already been trimmed')
        read=self.sequence
        if not linker:
            if self.space=="DNA": linkseq = P2_seq
            elif self.space == "CS": linkseq = P2_CS_seq[1:]
            linker = linker_oligos(linkseq)
    
        linker_len = len(linkseq)
    
        ##from max. possible overlap, check and take best
        max_ol = min([len(read), linker_len])
        for n in range(max_ol, 0, -1):
            if read[-n:] in linker:
                self.sequence = read[:-n]
                self.qual = self.qual[:-n]
                self.trimmed=True
                break
        return #self.sequence
    
    def trim_by_qual(self,phredCutoff=10):
        """iterative trimming of 3' end by quality cutoff (default = 10)"""
        bases = 0
        for q in self.qual[::-1]:
            if q>=phredCutoff: break
            else: bases+=1
        if bases >0:
            self.sequence = self.sequence[:-bases]
            self.qual = self.qual[:-bases]
        return
    
    def nuIDName(self):
        if self.space == "CS":
            tempString = CS2DNA(self.sequence)
        else:
            tempString = self.sequence
        nuID = misc.seq2nuID(tempString)
        self.name=nuID
        return
########################################################################
#Basic Iterators for SOLiD Data
########################################################################        
def CSFastaIterator(handle, matches=False):
    """
    Generator function to iterate over csfasta records in <handle>:
    Use in a loop to apply to each CSSeq record contained in a .csfasta file
    Input: record handle as obtained by handle = open(<file>,'r')
    Returns an iterator across CSSeq objects
    """
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break
    #Begin walk through csfasta records
    while True:
        if line[0] <>">":
            raise ValueError("Records in csfasta files should start with a '>' character")
        name = line[1:].rstrip()
        #if matches:
        parsedList = name.split(',')
        name = parsedList[0]
        matchList = parsedList[1:]
            #count = len(matchList)
        
        lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            lines.append(line.rstrip().replace(" ",""))
            line = handle.readline()
        
        #print matchList
        #Return record then continue
        newSeq = CSSeq(name,"".join(lines))
        if matches:
            newSeq.matches = matchList
        #if count != 0:
            #newSeq.count = count 
        yield newSeq
        
        if not line : return #StopIteration
    assert False, "Should not reach this line"

def QualIterator(handle):
    """Simple iterator over .qual files (try to use CompIter for paired analysis of .csfasta and .qual files)"""
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break
    while True:
        if line[0] <>">":
            raise ValueError("Records in .qual files should start with a '>' character")
        qual={}
        qual['name'] = line[1:].rstrip()
        #scores = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            try: 
                qual['scores']=map(int,line.rstrip().split())
            except ValueError:
                assert ValueError(" ".join([str(x) for x in qual['scores']]))
            line = handle.readline()
        
        yield qual
        
        if not line : return #StopIteration
    assert False, "Should not reach this line"

def CompIter(csfile,qualfile):
    """Takes a matched pair of .csfasta and .qual files and returns an iterator of CSSeq objects with a 'qual'
    attribute corresponding to the Phred scores for the particular bead
    """
    qualhandle=open(qualfile,'r')
    csfastahandle=open(csfile,'r')

    csiter=CSFastaIterator(csfastahandle)
    qualiter=QualIterator(qualhandle)

    for i in csiter:
        q=qualiter.next()   
        if q['name']==i.name:
            i.qual=q['scores']
            yield i
        else:
            assert ValueError ("It appears that the sequences don't match...have you modified the .csfasta or .qual files?")

def uniqueTableIterator(handle,trim=True):
    for line in handle:
        tokens = line.rstrip().split("\t")
        seq = CSSeq(tokens[0],tokens[0],readcount=int(tokens[1]))
        seq.nuIDName()
        if trim:
            seq.strip_solid_linker()
        yield seq

########################################################################
#MAQ-related functions
########################################################################
def SangerQualString(quals):
    """Takes a list of quality scores and returns the Sanger fastq string of quality values chr(val+33)"""
    return """""".join(chr(qual+33) for qual in quals)

def maqSeqTrim(CSsequence):
    """Converts a cs sequence to a maq-formatted string for alignment.
    Use only as intermediate to go into makeFastq.
    """
    return misc.mreplace(CSsequence[2:],chararray=['0','1','2','3'],newarray=['A','C','G','T'])

def makeFastq(csfile,qualfile,shortname,outdir="",split=-1,trim=False):
    """Takes a matched pair of SOLiD files (one .csfasta and one .qual) and prints a 'sanger formatted' fastq file
    to outfile
    Does not yet work for paired end reads.
    This is tweaked to make Bowtie-compatible fastq files...
    Use split to place only 'split' number of reads into each output file (as name implies).
    """
    iter = CompIter(csfile,qualfile)
    group = 1
    
    #Test to see if output directory is accessible and if not, it creates it. (This could be more streamlined)
    if outdir != "" and os.access(outdir, os.F_OK) is False:
        os.mkdir(outdir)
    outhand = open(outdir+shortname+"_"+str(group)+".fastq",'w')
    counter = 0
    for i in iter:
        if "." in i.sequence: continue
        counter += 1
        if trim:
            i.strip_solid_linker()
        print >>outhand, """@%s:%s/1\n%s\n+\n%s""" % (shortname,i.name[:-3],i.sequence,SangerQualString(i.qual))
        if split > 0 and counter%split == 0:
            group +=1
            outhand.close()
            outhand = open(outdir+shortname+"_"+str(group)+".fastq",'w')
    outhand.close()
    return

########################################################################
#Misc. Tools for SOLiD manipulation
########################################################################

def csfasta2fasta(fname):
    handle=open(fname,'r')
    iter=CSFastaIterator(handle)
    for i in iter:
        i.CS2DNA
        i.printFasta


def uniqueTable(dir=os.getcwd()):
    """
    Takes a all .csfasta files in a directory  and creates a dataframe (tab-delimited) 
    of counts for unique reads in each file
    REQUIRES A LARGE AMOUNT OF MEMORY AS READ INDEX IS STORED AS DICTIONARY!!!
    
    OR
    
    use awk script ~/bin/makeUnique.awk
    
    """
    files = os.listdir(dir)
    dict={}
    samples=[]
    for file in files:
        if file.rfind('.csfasta')>0:
            handle = open(file,'r')
            sample = file.rstrip('.csfasta')
            samples.append(sample)
            sys.stderr.write("Processing Sample %s...\n" % sample)
            iter = CSFastaIterator(handle)
            count = 0
            for i in iter:
                count = count + 1
                if not dict.has_key(i.sequence):
                    dict[i.sequence] = {}
                if not dict[i.sequence].has_key(sample):
                    dict[i.sequence][sample] = 1
                else:
                    dict[i.sequence][sample] = dict[i.sequence][sample] + 1
                if count % 100000 == 0:
                    sys.stderr.write("%d\n" % count)
            handle.close()
    sys.stderr.write("Retrieving Keys...\n")
    keys = dict.keys()
    sys.stderr.write("Sorting Keys...\n")
    keys.sort()
    sys.stderr.write("Writing to output...\n")
    samples.sort()
    print "#Sequence\t",
    print "\t".join(samples)
    for key in keys:
        print "%s\t" % key,
        #print dict[key]
        
        for sample in samples:
            if dict[key].has_key(sample):
                continue
            else:
                dict[key][sample] = 0
            
        #print dict[key]
        for sample in samples:
            print "%d\t" % dict[key][sample],
        print ""
        
def filterUnique(uniqueFile,minObs=5):
    """
    At this point, this function is specific to the H1U and H1NSC samples
    I need to change that
    """
    handle = open(uniqueFile,'r')
    count = 0
    Ufile = open('H1U.csfasta','w')
    NSCfile = open('H1NSC.csfasta','w')
    while True:
        line = handle.readline()
        count = count + 1
        if not line: break
        #if count == 10: break
        if line.startswith("#"):
            samples = line[1:].split()
            samples = samples[1:]
            #print samples
            continue
        values = line.split()
        readSeq = values.pop(0)
        values = map(int,values)
        ###
        #This is currently specific to H1U and H1NSC samples
        ###
        #total = sum(values)
        U = sum(values[0:2])
        NSC = sum(values[2:])
        if U>=minObs:
            Ufile.write(">%s_x%d\n%s\n" % (readSeq,U,readSeq))
        if NSC>=minObs:
            NSCfile.write(">%s_x%d\n%s\n" % (readSeq,NSC,readSeq))
    Ufile.close()
    NSCfile.close()
    
def CS2DNA(sequence):
    """
    Takes a colorspace sequence and converts it to DNA space
    """
    mapping = {}
    mapping["0"] = {"T":"T","A":"A","C":"C","G":"G"}
    mapping["1"] = {"T":"G","A":"C","C":"A","G":"T"}
    mapping["2"] = {"T":"C","A":"G","C":"T","G":"A"}
    mapping["3"] = {"T":"A","A":"T","C":"G","G":"C"}
    
    res = ''
    letter = ''
    
    for i in sequence:
        if (letter == ''):
            letter = res = i
            continue
        else:
            letter = mapping[i][letter]
        res = res + letter
    return res[1:]
