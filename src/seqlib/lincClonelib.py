#!/usr/bin/env python
'''
Created on Aug 19, 2010

Requirements:
    - primer3_core

@author: Loyal Goff

TODO:
- Add bed file output for primers as option
- Integrate a few more primer3 options into commandline
    * number of primers
    * GC adjustment
    * etc...
'''

#from Bio.Emboss import Primer3
from RNASeq import sequencelib,primer3lib
import subprocess,sys,getopt,os

help_message = '''
usage:
python lincClonelib.py [options] <fastaFile.fa>

options:
    -h or --help      Prints this helpful help message
    -o or --output    output file for pretty results (default = <fastaFile_primers.txt>
    -g                Add attB sites for gateway cloning
    -k                Keep tmp files
    -v                Verbose output
    -t                tab-delimited output (more machine readable)
'''

wiggleRoom = 40
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=36
clonePrimerSteps = [0,5,10,20,40,50]
attF = "GGGGACAAGTTTGTACAAAAAAGCAGGCT" #Sequence to be added to the forward primer for Gateway (TM) cloning
attR = "GGGGACCACTTTGTACAAGAAAGCTGGGT" #Sequence to be added to the reverse primer for Gateway (TM) cloning


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def runPrimer3(fastaFile,p3CloneSetFile="/n/rinn_data1/users/lgoff/utils/primer_design/P3_cloning_primer_settings.p3",p3PCRSetFile="/n/rinn_data1/users/lgoff/utils/primer_design/P3_qPCR_primer_settings.p3",p3InsituSetFile="/n/rinn_data1/users/lgoff/utils/primer_design/P3_insitu_probe_settings.p3",verbose=False,keepTmp=False):
    baseName = fastaFile.rstrip(".fa")
    iter = sequencelib.FastaIterator(open(fastaFile,'r'))
    cloneTmpFname = baseName+"_clone.p3in"
    cloneTmpHandle = open(cloneTmpFname,'w')
    qPCRTmpFname = baseName+"_qPCR.p3in"
    qPCRTmpHandle = open(qPCRTmpFname,'w')
    insituTmpFname = baseName+"_insitu.p3in"
    insituTmpHandle = open(insituTmpFname,'w')
    
    #Make Boulder-IO format...
    for i in iter:
        seqLength=len(i['sequence'])
        if seqLength-clonePrimerSteps[-1]<=PRIMER_MAX_SIZE:
            sys.stderr.write("%s sequence to short\n" % (i['name']))
            continue
        print >>qPCRTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n=" % (i['name'],i['sequence'])
        #print >>cloneTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nSEQUENCE_INCLUDED_REGION=1,%d\n=" % (i['name'],i['sequence'],len(i['sequence']))
        #print >>cloneTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,%d,%d,%d\n=" % (i['name'],i['sequence'],wiggleRoom,len(i['sequence'])-wiggleRoom,wiggleRoom)
        #print >>cloneTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nPRIMER_PRODUCT_SIZE_RANGE=%d-%d %d-%d %d-%d %d-%d %d-%d %d-%d\n=" % (i['name'],i['sequence'],len(i['sequence']),len(i['sequence']),len(i['sequence'])-5,len(i['sequence']),len(i['sequence'])-10,len(i['sequence']),len(i['sequence'])-20,len(i['sequence']),len(i['sequence'])-40,len(i['sequence']),len(i['sequence'])-50,len(i['sequence']))
        print >>cloneTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nSEQUENCE_INCLUDED_REGION=%d,%d\n=" % (i['name'],i['sequence'],1,len(i['sequence']))
        print >>insituTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n=" % (i['name'],i['sequence'])        
        
    qPCRTmpHandle.close()
    cloneTmpHandle.close()
    insituTmpHandle.close()
    
    P3Command = "primer3_core -p3_settings_file=%s -output=%s.p3out %s"
    #P3Command = "primer3_core -format_output -p3_settings_file=%s -output=%s.p3out %s"
    
    if verbose:
        sys.stderr.write("Designing qPCR Primers...\n")
    qpcr = subprocess.Popen(P3Command % (p3PCRSetFile,baseName+"_qPCR",qPCRTmpFname),shell=True)
    if verbose:
        sys.stderr.write("Designing Cloning Primers...\n")
    cloning = subprocess.Popen(P3Command % (p3CloneSetFile,baseName+"_cloning",cloneTmpFname),shell=True)
    if verbose:
        sys.stderr.write("Designing InSitu Primers...\n")
    insitu = subprocess.Popen(P3Command % (p3InsituSetFile,baseName+"_insitu",insituTmpFname),shell=True)
    qpcr.wait()
    cloning.wait()
    insitu.wait()
    if not keepTmp:
        os.remove(cloneTmpFname)
        os.remove(qPCRTmpFname)
        os.remove(insituTmpFname)
    return (baseName+"_qPCR.p3out",baseName+"_cloning.p3out",baseName+"_insitu.p3out")
    
def test():
    fastaFile="lincSFPQ.fa"
    qPCR,cloning = runPrimer3(fastaFile)
    return

def parsePrimer3(p3OutFile):
    handle = open(p3OutFile,'r')
    iter = primer3lib.parse(handle)
    for record in iter:
        yield record

def printqPCR(p3outFile,outHandle):
    recordIter = parsePrimer3(p3outFile)
    print >>outHandle, "######################\n# qPCR Primers\n######################"
    for record in recordIter:
        print >>outHandle, "%s" % record.sequenceID
        if len(record.primers)<1:
            print >>outHandle, "\tNo acceptable qPCR primers were found."
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                
                fwdStr = "\t%d) Amplicon Size: %d\n\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc)
                revStr = "\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, fwdStr
                print >>outHandle, revStr
                print >>outHandle, ""
        print >>outHandle, "--------------------------------"

def printqPCRTabDelim(p3outFile,outHandle):
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "######################\n# qPCR Primers\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print >>outHandle, "%s\tqPCR\t%s" % (record.sequenceID,'No acceptable qPCR primers were found.')
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                outStr = "%s\tqPCR\t%d\t%d\t%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%d\t%d\t%0.2f\t%0.2f" % (record.sequenceID,primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc,revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, outStr


def printCloning(p3outFile,outHandle,gateway=False):
    recordIter = parsePrimer3(p3outFile)
    print >>outHandle, "\n######################\n# Cloning Primers\n######################"
    for record in recordIter:
        print >>outHandle, "%s" % record.sequenceID
        if len(record.primers)<1:
            print >>outHandle, "\tNo acceptable Cloning primers were found."
            continue
        else:
            for primer in record.primers:
                if gateway:
                    fwdSeq = attF+primer.forward_seq
                    revSeq = attR+primer.reverse_seq
                    gatewayStr = "Gateway"
                else:
                    fwdSeq = primer.forward_seq
                    revSeq = primer.reverse_seq
                    gatewayStr = ""
                fwdStr = "\t%d) Amplicon Size: %d\t%s\n\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (primer.number,primer.product_size,gatewayStr,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc)
                revStr = "\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, fwdStr
                print >>outHandle, revStr
                print >>outHandle, ""
        print >>outHandle, "--------------------------------"

def printCloningTabDelim(p3outFile,outHandle,gateway=False):
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "\n######################\n# Cloning Primers\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print >>outHandle, "%s\tCloning\t%s" % (record.sequenceID,'No acceptable primers were found.')
            continue
        else:
            for primer in record.primers:
                if gateway:
                    fwdSeq = attF+primer.forward_seq
                    revSeq = attR+primer.reverse_seq
                    gatewayStr = "Gateway"
                else:
                    fwdSeq = primer.forward_seq
                    revSeq = primer.reverse_seq
                    gatewayStr = ""
                outStr = "%s\tCloning\t%d\t%d\t%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%d\t%d\t%0.2f\t%0.2f" % (record.sequenceID,primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc,revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, outStr

def printInsitu(p3outFile,outHandle):
    recordIter = parsePrimer3(p3outFile)
    print >>outHandle, "######################\n# InSitu Primers\n######################"
    for record in recordIter:
        print >>outHandle, "%s" % record.sequenceID
        if len(record.primers)<1:
            print >>outHandle, "\tNo acceptable InSitu primers were found."
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                
                fwdStr = "\t%d) Amplicon Size: %d\n\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc)
                revStr = "\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, fwdStr
                print >>outHandle, revStr
                print >>outHandle, ""
        print >>outHandle, "--------------------------------"

def printInsituTabDelim(p3outFile,outHandle):
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "######################\n# qPCR Primers\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print >>outHandle, "%s\tInSitu\t%s" % (record.sequenceID,'No acceptable InSitu primers were found.')
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                outStr = "%s\tInSitu\t%d\t%d\t%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%d\t%d\t%0.2f\t%0.2f" % (record.sequenceID,primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc,revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, outStr

def printInsitu(p3outFile,outHandle):
    recordIter = parsePrimer3(p3outFile)
    print >>outHandle, "######################\n# InSitu Primers\n######################"
    for record in recordIter:
        print >>outHandle, "%s" % record.sequenceID
        if len(record.primers)<1:
            print >>outHandle, "\tNo acceptable InSitu primers were found."
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                
                fwdStr = "\t%d) Amplicon Size: %d\n\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc)
                revStr = "\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, fwdStr
                print >>outHandle, revStr
                print >>outHandle, ""
        print >>outHandle, "--------------------------------"

def printInsituTabDelim(p3outFile,outHandle):
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "######################\n# ASO Candidates\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print >>outHandle, "%s\tASO\t%s" % (record.sequenceID,'No acceptable ASO candidates were found.')
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                outStr = "%s\tInSitu\t%d\t%d\t%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%d\t%d\t%0.2f\t%0.2f" % (record.sequenceID,primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc,revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print >>outHandle, outStr

def main(argv=None): 
    if argv is None:
        argv = sys.argv
    task = 'qpcr'
    verbose = False
    outFile = None
    gateway = False
    keepTmp = False
    tabDelim = False
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hto:vgk", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option == "-g":
                gateway = True
            if option == "-k":
                keepTmp=True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                outFile = value
            if option == "-t":
                tabDelim = True
        try:
            assert len(args)==1
            fname=args[0]
        except:
            raise Usage(help_message)
        if outFile == None:
            outFile = fname.rstrip(".fa")+"_primers.txt"
        outHandle = open(outFile,'w')
        qPCR,cloning,insitu = runPrimer3(fname,verbose=verbose,keepTmp=keepTmp)
        if tabDelim:
            print >>outHandle, "sequenceID\tPrimer Type\tPrimer number\tProduct_size\tFwdSeq\tForward start\tLength Fwd\tFwd Tm\tFwd GC\tRevSeq\tRev start\tLength Rev\tRev Tm\tRev GC"
            printqPCRTabDelim(qPCR,outHandle)
            printCloningTabDelim(cloning,outHandle,gateway=gateway)
            printInsituTabDelim(insitu,outHandle)
        else:
            printqPCR(qPCR,outHandle)
            printCloning(cloning,outHandle,gateway=gateway)
            printInsitu(insitu,outHandle)
        if not keepTmp:
            os.remove(qPCR)
            os.remove(cloning)
            os.remove(insitu)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        sys.exit()
    

if __name__ == "__main__":
    sys.exit(main())
