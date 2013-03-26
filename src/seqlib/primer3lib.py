'''
Created on Sep 9, 2010

Handles primer3 running and parsing output

primer3 >= v2.2

@author: lgoff
'''
import sys,subprocess
from RNASeq import sequencelib


class Record(object):
    '''
    Represent information from a primer3 run finding primers.
    
    Members:
        - sequenceID = value of SEQUENCE_ID field from primer3 record
        - sequence = value of SEQUENCE_TEMPLATE field 
        - primers = list of Primer objects describing primer pairs for this target sequence.
        - comments = the comment line(s) for the record
        - attributes = other global parameters relevant to the record as a whole and not just a primer
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.sequenceID = ""
        self.sequence = ""
        self.comments = ""
        self.primers = []
        self.attributes = {}
    
    def __iter__(self):
        return iter(self.primers)
    
    def __repr__(self):
        return "%s: %d primer pair(s)" % (self.sequenceID,len(self.primers))
    
class Primer(object):
    '''
    A primer set designed by Primer3
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.sequenceID=""
        self.number = 0
        self.size = 0
        self.forward_seq = ''
        self.forward_start = ''
        self.forward_length = ''
        self.forward_tm = 0.0
        self.forward_gc = 0.0
        self.reverse_seq = ''
        self.reverse_start = 0
        self.reverse_length = 0
        self.reverse_tm = 0.0
        self.reverse_gc = 0.0
        self.product_size = 0
    
    def __repr__(self):
        return "%s_%d\n\tFwd: %s\tRev: %s" % (self.sequenceID,self.number,self.forward_seq, self.reverse_seq)
    
def parse(handle):
    recordLines = []
    while True:
        line = handle.readline().rstrip()
        if not line: raise StopIteration
        if not line == "=":
            recordLines.append(line)
            continue
        else:
            recordLines = [x.split("=") for x in recordLines]
            recordDict = dict(zip([x[0] for x in recordLines],[x[1] for x in recordLines]))
            rdKeys = recordDict.keys()
            record = Record()
            record.sequenceID = recordDict['SEQUENCE_ID']
            record.sequence = recordDict['SEQUENCE_TEMPLATE']
            try:
                nPrimers = int(recordDict['PRIMER_PAIR_NUM_RETURNED'])
            except KeyError:
                nPrimers=0
            for i in xrange(nPrimers):
                primer = Primer()
                primer.sequenceID = record.sequenceID
                primer.number = i+1
                primer.size = int(recordDict['PRIMER_PAIR_%d_PRODUCT_SIZE' % i])
                primer.forward_seq = recordDict['PRIMER_LEFT_%d_SEQUENCE' % i]
                primer.forward_start = int(recordDict['PRIMER_LEFT_%d' % i].split(",")[0])
                primer.forward_length = int(recordDict['PRIMER_LEFT_%d' % i].split(",")[1])
                primer.forward_tm = float(recordDict['PRIMER_LEFT_%d_TM' % i])
                primer.forward_gc = float(recordDict['PRIMER_LEFT_%d_GC_PERCENT' % i])
                primer.reverse_seq = recordDict['PRIMER_RIGHT_%d_SEQUENCE' % i]
                primer.reverse_start = int(recordDict['PRIMER_RIGHT_%d' % i].split(",")[0])
                primer.reverse_length = int(recordDict['PRIMER_RIGHT_%d' % i].split(",")[1])
                primer.reverse_tm = float(recordDict['PRIMER_RIGHT_%d_TM' % i])
                primer.reverse_gc = float(recordDict['PRIMER_RIGHT_%d_GC_PERCENT' % i])
                primer.product_size = int(recordDict['PRIMER_PAIR_%d_PRODUCT_SIZE' % i])
                record.primers.append(primer)
            yield record
            recordLines = []

#######
#Context specific runs
#######
def runPrimer3(fastaFile,task="qpcr",p3CloneSetFile="/seq/compbio-hp/lgoff/lincRNAs/primer_design/P3_cloning_primer_settings.p3",p3PCRSetFile="/seq/compbio-hp/lgoff/lincRNAs/primer_design/P3_qPCR_primer_settings.p3"):
    """Task can be either 'qpcr' or 'cloning'"""
    
    baseName = fastaFile.rstrip(".fa")
    iter = sequencelib.FastaIterator(open(fastaFile,'r'))
    tmpFname = baseName+".p3in"
    tmpHandle = open(tmpFname,'w')
    
    #Make Boulder-IO format...
    for i in iter:
        myString = "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n" % (i['name'],i['sequence'])
        if task == "cloning":
            myString += "SEQUENCE_INCLUDED_REGION=1,%d\n" % (i['name'],i['sequence'],len(i['sequence']))
        myString += "="
        print >>tmpHandle, myString
    tmpHandle.close()
    
    P3Command = "primer3_core -p3_settings_file=%s -output=%s.p3out %s"
    
    sys.stderr.write("Designing Primers...\n")
    if task == "qpcr":
        subprocess.Popen(P3Command % (p3PCRSetFile,baseName+"_qPCR",tmpFname),shell=True)
    elif task == "cloning":
        subprocess.Popen(P3Command % (p3CloneSetFile,baseName+"_cloning",tmpFname),shell=True)
    return baseName+".p3out"
