#!/usr/bin/env python
"""Primer design pipeline for lincRNA cloning, qPCR, and in situ hybridisation.

Wraps the primer3_core command-line tool to design three classes of primers
from FASTA sequences: cloning primers (with optional Gateway attB flanks),
qPCR primers, and in situ hybridisation probe primers.  Output can be
formatted as human-readable text or as tab-delimited tables for downstream
processing.

Requirements:
    - primer3_core executable on PATH

Usage::

    python lincClonelib.py [options] <fastaFile.fa>
"""

#from Bio.Emboss import Primer3
import getopt
import os
import subprocess
import sys

from RNASeq import primer3lib, sequencelib

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
    """Exception raised for command-line usage errors in lincClonelib.

    Attributes:
        msg: Human-readable explanation of the error or the help message.
    """
    def __init__(self, msg):
        """Initialises a Usage exception with an error message.

        Args:
            msg: Human-readable error or help text.
        """
        self.msg = msg

def runPrimer3(fastaFile,p3CloneSetFile="/n/rinn_data1/users/lgoff/utils/primer_design/P3_cloning_primer_settings.p3",p3PCRSetFile="/n/rinn_data1/users/lgoff/utils/primer_design/P3_qPCR_primer_settings.p3",p3InsituSetFile="/n/rinn_data1/users/lgoff/utils/primer_design/P3_insitu_probe_settings.p3",verbose=False,keepTmp=False):
    """Runs primer3_core to design qPCR, cloning, and in situ primers from a FASTA file.

    Creates three Boulder-IO input files from the FASTA sequences and launches
    three parallel primer3_core processes (one per primer type), each with its
    own settings file.  Waits for all processes to complete before returning.

    Args:
        fastaFile: Path to a FASTA file of sequences to design primers for.
            Sequences shorter than clonePrimerSteps[-1] + PRIMER_MAX_SIZE
            bases are skipped for cloning design.
        p3CloneSetFile: Path to a primer3 settings file for cloning primers.
        p3PCRSetFile: Path to a primer3 settings file for qPCR primers.
        p3InsituSetFile: Path to a primer3 settings file for in situ primers.
        verbose: If True, write progress messages to stderr (default: False).
        keepTmp: If True, retain the temporary Boulder-IO input files after
            the run (default: False).

    Returns:
        A tuple of three strings: (qPCR_output_path, cloning_output_path,
        insitu_output_path) giving the paths to the primer3 output files.
    """
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
        print("SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n=" % (i['name'],i['sequence']), file=qPCRTmpHandle)
        #print >>cloneTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nSEQUENCE_INCLUDED_REGION=1,%d\n=" % (i['name'],i['sequence'],len(i['sequence']))
        #print >>cloneTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,%d,%d,%d\n=" % (i['name'],i['sequence'],wiggleRoom,len(i['sequence'])-wiggleRoom,wiggleRoom)
        #print >>cloneTmpHandle, "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nPRIMER_PRODUCT_SIZE_RANGE=%d-%d %d-%d %d-%d %d-%d %d-%d %d-%d\n=" % (i['name'],i['sequence'],len(i['sequence']),len(i['sequence']),len(i['sequence'])-5,len(i['sequence']),len(i['sequence'])-10,len(i['sequence']),len(i['sequence'])-20,len(i['sequence']),len(i['sequence'])-40,len(i['sequence']),len(i['sequence'])-50,len(i['sequence']))
        print("SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nSEQUENCE_INCLUDED_REGION=%d,%d\n=" % (i['name'],i['sequence'],1,len(i['sequence'])), file=cloneTmpHandle)
        print("SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n=" % (i['name'],i['sequence']), file=insituTmpHandle)

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
    """Smoke test for runPrimer3 using a hard-coded FASTA file.

    Calls runPrimer3 on 'lincSFPQ.fa' and returns nothing.  Intended for
    interactive testing only.
    """
    fastaFile="lincSFPQ.fa"
    qPCR,cloning = runPrimer3(fastaFile)
    return

def parsePrimer3(p3OutFile):
    """Yields parsed primer3 Record objects from a primer3 output file.

    Opens the specified output file and delegates parsing to primer3lib.parse,
    yielding one Record object per sequence entry.

    Args:
        p3OutFile: Path to a primer3 output file (Boulder-IO format).

    Yields:
        primer3lib.Record objects, each containing the sequenceID, template
        sequence, and a list of Primer objects.
    """
    handle = open(p3OutFile,'r')
    iter = primer3lib.parse(handle)
    for record in iter:
        yield record

def printqPCR(p3outFile,outHandle):
    """Writes qPCR primer results in human-readable format.

    Parses primer3 output and writes a formatted, multi-line report of qPCR
    primer pairs grouped by sequence ID.  If no acceptable primers were found
    for a sequence, a placeholder message is printed.

    Args:
        p3outFile: Path to a primer3 qPCR output file.
        outHandle: Writable file-like object to receive the formatted output.
    """
    recordIter = parsePrimer3(p3outFile)
    print("######################\n# qPCR Primers\n######################", file=outHandle)
    for record in recordIter:
        print("%s" % record.sequenceID, file=outHandle)
        if len(record.primers)<1:
            print("\tNo acceptable qPCR primers were found.", file=outHandle)
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq

                fwdStr = "\t%d) Amplicon Size: %d\n\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc)
                revStr = "\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print(fwdStr, file=outHandle)
                print(revStr, file=outHandle)
                print("", file=outHandle)
        print("--------------------------------", file=outHandle)

def printqPCRTabDelim(p3outFile,outHandle):
    """Writes qPCR primer results in tab-delimited format.

    Parses primer3 output and writes one line per primer pair with columns:
    sequenceID, primer type ('qPCR'), primer number, product size, forward
    sequence, forward start, forward length, forward Tm, forward GC, reverse
    sequence, reverse start, reverse length, reverse Tm, reverse GC.

    Args:
        p3outFile: Path to a primer3 qPCR output file.
        outHandle: Writable file-like object to receive the tab-delimited
            output.
    """
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "######################\n# qPCR Primers\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print("%s\tqPCR\t%s" % (record.sequenceID,'No acceptable qPCR primers were found.'), file=outHandle)
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                outStr = "%s\tqPCR\t%d\t%d\t%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%d\t%d\t%0.2f\t%0.2f" % (record.sequenceID,primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc,revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print(outStr, file=outHandle)


def printCloning(p3outFile,outHandle,gateway=False):
    """Writes cloning primer results in human-readable format.

    Parses primer3 output and writes a formatted, multi-line report of
    cloning primer pairs grouped by sequence ID.  When gateway is True,
    Gateway attB sequences are prepended to the forward and reverse primers
    and 'Gateway' is noted in the output.

    Args:
        p3outFile: Path to a primer3 cloning output file.
        outHandle: Writable file-like object to receive the formatted output.
        gateway: If True, prepend attF to forward and attR to reverse primers
            for Gateway cloning (default: False).
    """
    recordIter = parsePrimer3(p3outFile)
    print("\n######################\n# Cloning Primers\n######################", file=outHandle)
    for record in recordIter:
        print("%s" % record.sequenceID, file=outHandle)
        if len(record.primers)<1:
            print("\tNo acceptable Cloning primers were found.", file=outHandle)
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
                print(fwdStr, file=outHandle)
                print(revStr, file=outHandle)
                print("", file=outHandle)
        print("--------------------------------", file=outHandle)

def printCloningTabDelim(p3outFile,outHandle,gateway=False):
    """Writes cloning primer results in tab-delimited format.

    Parses primer3 output and writes one line per primer pair with columns:
    sequenceID, primer type ('Cloning'), primer number, product size, forward
    sequence, forward start, forward length, forward Tm, forward GC, reverse
    sequence, reverse start, reverse length, reverse Tm, reverse GC.  When
    gateway is True, attB sequences are prepended to the primer sequences.

    Args:
        p3outFile: Path to a primer3 cloning output file.
        outHandle: Writable file-like object to receive the tab-delimited
            output.
        gateway: If True, prepend attF to forward and attR to reverse primers
            (default: False).
    """
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "\n######################\n# Cloning Primers\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print("%s\tCloning\t%s" % (record.sequenceID,'No acceptable primers were found.'), file=outHandle)
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
                print(outStr, file=outHandle)

def printInsitu(p3outFile,outHandle):
    """Writes in situ hybridisation primer results in human-readable format.

    Parses primer3 output and writes a formatted, multi-line report of in situ
    probe primer pairs grouped by sequence ID.

    Args:
        p3outFile: Path to a primer3 in situ output file.
        outHandle: Writable file-like object to receive the formatted output.
    """
    recordIter = parsePrimer3(p3outFile)
    print("######################\n# InSitu Primers\n######################", file=outHandle)
    for record in recordIter:
        print("%s" % record.sequenceID, file=outHandle)
        if len(record.primers)<1:
            print("\tNo acceptable InSitu primers were found.", file=outHandle)
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq

                fwdStr = "\t%d) Amplicon Size: %d\n\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc)
                revStr = "\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print(fwdStr, file=outHandle)
                print(revStr, file=outHandle)
                print("", file=outHandle)
        print("--------------------------------", file=outHandle)

def printInsituTabDelim(p3outFile,outHandle):
    """Writes in situ hybridisation primer results in tab-delimited format.

    Parses primer3 output and writes one line per primer pair with columns:
    sequenceID, primer type ('InSitu'), primer number, product size, forward
    sequence, forward start, forward length, forward Tm, forward GC, reverse
    sequence, reverse start, reverse length, reverse Tm, reverse GC.

    Args:
        p3outFile: Path to a primer3 in situ output file.
        outHandle: Writable file-like object to receive the tab-delimited
            output.
    """
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "######################\n# qPCR Primers\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print("%s\tInSitu\t%s" % (record.sequenceID,'No acceptable InSitu primers were found.'), file=outHandle)
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                outStr = "%s\tInSitu\t%d\t%d\t%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%d\t%d\t%0.2f\t%0.2f" % (record.sequenceID,primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc,revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print(outStr, file=outHandle)

def printInsitu(p3outFile,outHandle):
    """Writes in situ hybridisation primer results in human-readable format (second definition).

    Duplicate of the earlier printInsitu definition; this version is the one
    that Python will actually use at runtime.  Parses primer3 output and writes
    a formatted, multi-line report of in situ probe primer pairs grouped by
    sequence ID.

    Args:
        p3outFile: Path to a primer3 in situ output file.
        outHandle: Writable file-like object to receive the formatted output.
    """
    recordIter = parsePrimer3(p3outFile)
    print("######################\n# InSitu Primers\n######################", file=outHandle)
    for record in recordIter:
        print("%s" % record.sequenceID, file=outHandle)
        if len(record.primers)<1:
            print("\tNo acceptable InSitu primers were found.", file=outHandle)
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq

                fwdStr = "\t%d) Amplicon Size: %d\n\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc)
                revStr = "\t\t%s\tStart: %d\tLength: %d\tTm: %0.2f\tGC: %0.2f" % (revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print(fwdStr, file=outHandle)
                print(revStr, file=outHandle)
                print("", file=outHandle)
        print("--------------------------------", file=outHandle)

def printInsituTabDelim(p3outFile,outHandle):
    """Writes ASO / in situ primer results in tab-delimited format (second definition).

    Duplicate of the earlier printInsituTabDelim definition; this version
    overrides the first at runtime.  Parses primer3 output for in situ /
    ASO candidates and writes one tab-delimited line per primer pair with
    an 'InSitu' type column.  When no candidates are found, writes an 'ASO'
    type placeholder line.

    Args:
        p3outFile: Path to a primer3 output file.
        outHandle: Writable file-like object to receive the tab-delimited
            output.
    """
    recordIter = parsePrimer3(p3outFile)
    #print >>outHandle, "######################\n# ASO Candidates\n######################"
    for record in recordIter:
        if len(record.primers)<1:
            print("%s\tASO\t%s" % (record.sequenceID,'No acceptable ASO candidates were found.'), file=outHandle)
            continue
        else:
            for primer in record.primers:
                #This is in place to extend the primer sequences with Restriction Sites at a later date if necessary...
                fwdSeq = primer.forward_seq
                revSeq = primer.reverse_seq
                outStr = "%s\tInSitu\t%d\t%d\t%s\t%d\t%d\t%0.2f\t%0.2f\t%s\t%d\t%d\t%0.2f\t%0.2f" % (record.sequenceID,primer.number,primer.product_size,fwdSeq,primer.forward_start,len(fwdSeq),primer.forward_tm,primer.forward_gc,revSeq,primer.reverse_start,len(revSeq),primer.reverse_tm,primer.reverse_gc)
                print(outStr, file=outHandle)

def main(argv=None):
    """Command-line entry point for the lincRNA primer design pipeline.

    Parses command-line options, runs primer3 via runPrimer3, and writes
    formatted primer output (human-readable or tab-delimited) to the output
    file.  Cleans up temporary primer3 output files unless keepTmp is set.

    Args:
        argv: List of command-line argument strings.  Defaults to sys.argv
            when None.

    Raises:
        SystemExit: On usage errors or when --help is requested.
    """
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
        except getopt.error as msg:
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
            print("sequenceID\tPrimer Type\tPrimer number\tProduct_size\tFwdSeq\tForward start\tLength Fwd\tFwd Tm\tFwd GC\tRevSeq\tRev start\tLength Rev\tRev Tm\tRev GC", file=outHandle)
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

    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg), file=sys.stderr)
        print("\t for help use --help", file=sys.stderr)
        sys.exit()


if __name__ == "__main__":
    sys.exit(main())
