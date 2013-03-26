'''
Created on Sep 19, 2012

Script to create gibson assembly fragments for ordering from a fasta file.

@author: lgoff
'''
#Imports
from RNASeq import sequencelib
from RNASeq.misc import pp
import getopt,sys,os


#Fixed attributes
attF = "GGGGACAAGTTTGTACAAAAAAGCAGGCT" #Sequence to be added to the forward primer for Gateway (TM) cloning
attR = "GGGGACCACTTTGTACAAGAAAGCTGGGT" #Sequence to be added to the reverse primer for Gateway (TM) cloning

#Error trapping
help_message = '''
usage:
python gibson.py [options] <fastaFile.fa>

options:
    -h or --help      Prints this helpful help message
    -o or --output    output file for pretty results (default = <fastaFile_primers.txt>
    -g                Add attB sites for gateway cloning
    -f                Fragment size (default: 500bp)
    -v                Verbose output
    -s                overhang size (default: 20bp)
    -t                tab-delimited output (more machine readable)
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def gibson(fname,gateway=True,fragSize=500,overhangSize=20):
    res = {}
    
    #Fasta file handle
    handle = open(fname,'r')
    iter = sequencelib.FastaIterator(handle)
    
    #Iterate over records in input fasta file
    for i in iter:
        fragments = []
        seq = i['sequence'].upper()
        if gateway:
            seq = attF + seq + sequencelib.rcomp(attR)
        curpos = 0
        length = int(len(seq)-1)
        while curpos < length:
            if curpos < 0:
                curpos = 0
            fragStart = curpos
            fragEnd = min(curpos+fragSize,length)
            #print "%d\t%d" % (fragStart,fragEnd)
            fragSeq = seq[int(fragStart):int(fragEnd)]
            fragments.append(fragSeq)
            curpos = curpos+fragSize-overhangSize
        res[i['name']]=fragments
    
    return res

def printGibson(fragDict,outHandle):
    for k in fragDict.keys():
        print >>outHandle, "%s:" % k
        blockCount = 0
        for fragment in fragDict[k]:
            blockCount += 1
            print >>outHandle,"%s_block%d\t%s" % (k,blockCount,fragment)
        print >>outHandle, "\n"
    
    

##############
# Main
##############
def main(argv=None):
    if argv is None:
        argv = sys.argv
    verbose = False
    outFile = None
    gateway = False
    keepTmp = False
    tabDelim = False
    overhangSize = 20
    fragSize = 500
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hto:vs:gf:k", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option == "-g":
                gateway = True
            if option == "-f":
                fragSize == value
            if option == "-k":
                keepTmp=True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                outFile = value
            if option == "-s":
                overhangSize=value
            if option == "-t":
                tabDelim = True
        try:
            assert len(args)==1
            fname=args[0]
        except:
            raise Usage(help_message)
        if outFile == None:
            outFile = fname.rstrip(".fa")+"_gibson.txt"
        outHandle = open(outFile,'w')
        
        #Put actual function call here...
        fragDict = gibson(fname,gateway=gateway,fragSize=fragSize,overhangSize=overhangSize)
        #pp(fragDict)
        printGibson(fragDict,outHandle)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        sys.exit()

if __name__ == "__main__":
    sys.exit(main())