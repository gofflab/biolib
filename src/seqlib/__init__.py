#!/usr/bin/env python
"""
Implementation of my short RNA Sequencing pipeline:
    Currently only for SHRiMP
    
    Usage: RNASeq.py -i input_file.csfasta -s shrimp_dir -o analysis_dir -a shrimp
    
    TODO:  
        -Adapt for MAQ and/or BOWTIE
        -Add module(s) for whole transcriptome analysis
            -exons
            -gene intersections
"""
#from shrimp import *
import sys,os,glob,getopt


def usage():
    pass

def main():
    try:
        opts,args = getopt.getopt(sys.argv[1:],'hvi:o:s:n:a',['help','verbose'])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    verbose = False
    aligner = 'shrimp'
    shrimpdir = os.getcwd()
    analyisdir = os.getcwd()
    samplename = "misc"
    
    for o,a in opts:
        if o == '-v':
            verbose = True
        elif o in ('-h','--help'):
            usage()
            sys.exit()
        elif o == '-i':
            fname = a
        elif o == '-s':
            shrimpdir = a
        elif o == '-o':
            analysisdir = a
        elif o == '-n':
            samplename = a
        elif o == 'a':
            aligner = a
        else:
            assert False, "Unhandled option"
    #Option checking
    if not fname.endswith('.csfasta'):
        print "Input file must be .csfasta format (appropriate extension required)"
        sys.exit(2)
    
    #Make directory structure for project    
    os.makedirs(shrimpdir+"/reads")
    os.makedirs(shrimpdir+"/results/split")
    if not analysisdir == os.getcwd():
        os.makedirs(analysisdir)
    
    #Split input .csfasta file
    sys.stderr.write("Splitting input file into reads directory")
    split_shrimp(fname,shrimpdir,binSize=1000)
    
    #TODO what the hell do I do with the LSF jobs after submission?
    

if __name__=="__main__":
    main()
    
    