'''
Created on Dec 15, 2009

Python tools for bowtie in colorspace (on Broad cluster)

@author: lgoff

Pipeline:
1) make cs .fastq file from .csfasta and .qual (see solid.py or use makeFastq.py in scripts)
2) Align reads in .fastq file using bowtie and specify SAM output with -S flag
3) Enjoy your alignments!


Example commandline:
bowtie -C -t -S -n 2 -k 1 -p 4 --best /seq/compbio-hp/lgoff/genomes/hg18/bowtie/hg18_c head0073_20090130_1Uppsala1_Upp1_F3_no_header.csfasta >head0073_20090130_1Uppsala1_Upp1_F3_bowtie.sam 2>bowtie.err
"""

'''
############
#Imports
############
import solid
import sys,os
############
#Constants
############

hg18_bowtieIndex = '/seq/compbio-hp/lgoff/genomes/hg18/bowtie/hg18_c'

########
#Prep CS for Bowtie
########

def prepBowtie(csfile,qualfile,shortname,basedir,split=100000,readsdir="fastq/",resultsdir="results/"):
    if not csfile.endswith('.csfasta'):
        raise ValueError("prepBowtie requires a .csfasta file")
    if not qualfile.endswith('.qual'):
        raise ValueError("prepBowtie requires a .qual file")
    #Make .fastq files
    sys.stderr.write("Making .fastq files...\n")
    solid.makeFastq(csfile,qualfile,shortname,outdir=readsdir,split=split)
    
    #Make resultsdir
    if os.access(resultsdir, os.F_OK) is False:
        os.mkdir(resultsdir)
    return

def runBowtie(queue="broad",cwd=os.getcwd(),outDir = "../results/"):
    files = os.listdir(cwd)
    for file in files:
        if file.endswith(".fastq"):
            basename = file.rstrip(".fastq")
            call = """bsub -q %s -P compbiofolk -o /dev/null -N "bowtie -C -t -S -n 2 -k 1 --best %s %s >%s%s.sam 2>%s%s.err" """ % (queue, hg18_bowtieIndex,file, outDir, basename, outDir, basename)  
            os.system(call)
            
    
    