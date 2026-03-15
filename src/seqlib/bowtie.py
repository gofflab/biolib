'''
Python tools for running Bowtie in colorspace mode on the Broad Institute cluster.

Provides helpers for preparing SOLiD colorspace reads for Bowtie alignment
and for submitting alignment jobs to an LSF cluster.  The pipeline is:

1. Make colorspace FASTQ files from ``.csfasta`` and ``.qual`` files
   (see ``solid.py`` or ``makeFastq.py``).
2. Align reads with ``bowtie`` using the ``-C`` (colorspace) and ``-S``
   (SAM output) flags.
3. Process the resulting SAM/BAM files with the tools in ``mySam.py`` or
   ``bwa.py``.

The module-level constant ``hg18_bowtieIndex`` points to the colorspace
Bowtie index used by the original author; update it for other references.

Created on Dec 15, 2009

@author: lgoff

Example commandline:
bowtie -C -t -S -n 2 -k 1 -p 4 --best /seq/compbio-hp/lgoff/genomes/hg18/bowtie/hg18_c head0073_20090130_1Uppsala1_Upp1_F3_no_header.csfasta >head0073_20090130_1Uppsala1_Upp1_F3_bowtie.sam 2>bowtie.err
'''
############
#Imports
############
import os
import sys

from . import solid

############
#Constants
############

hg18_bowtieIndex = '/seq/compbio-hp/lgoff/genomes/hg18/bowtie/hg18_c'

########
#Prep CS for Bowtie
########

def prepBowtie(csfile,qualfile,shortname,basedir,split=100000,readsdir="fastq/",resultsdir="results/"):
    """Prepare SOLiD colorspace reads for a Bowtie alignment run.

    Validates input file extensions, generates split FASTQ files from the
    colorspace FASTA and quality files using ``solid.makeFastq``, and
    creates the results output directory if it does not already exist.

    Args:
        csfile: Path to the SOLiD colorspace FASTA file (must end with
            ``.csfasta``).
        qualfile: Path to the quality score file (must end with ``.qual``).
        shortname: Base name used when naming the output FASTQ files.
        basedir: Base directory for the project (currently unused in the
            function body but reserved for future use).
        split: Maximum number of reads per split FASTQ file.  Defaults to
            100000.
        readsdir: Subdirectory path (relative to cwd) into which the FASTQ
            files are written.  Defaults to ``'fastq/'``.
        resultsdir: Subdirectory path (relative to cwd) that will receive
            Bowtie output.  Created if absent.  Defaults to ``'results/'``.

    Raises:
        ValueError: If ``csfile`` does not end with ``.csfasta``.
        ValueError: If ``qualfile`` does not end with ``.qual``.
    """
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
    """Submit colorspace Bowtie alignment jobs to an LSF cluster.

    Scans ``cwd`` for files ending in ``.fastq`` and submits one LSF
    ``bsub`` job per file.  Each job runs Bowtie in colorspace mode
    (``-C``), reporting a single best-alignment SAM file per input.

    Args:
        queue: LSF queue name to submit jobs to.  Defaults to ``'broad'``.
        cwd: Directory to scan for ``.fastq`` files.  Defaults to the
            current working directory at import time.
        outDir: Directory (relative or absolute) into which the SAM and
            error files are written.  Defaults to ``'../results/'``.
    """
    files = os.listdir(cwd)
    for file in files:
        if file.endswith(".fastq"):
            basename = file.rstrip(".fastq")
            call = """bsub -q %s -P compbiofolk -o /dev/null -N "bowtie -C -t -S -n 2 -k 1 --best %s %s >%s%s.sam 2>%s%s.err" """ % (queue, hg18_bowtieIndex,file, outDir, basename, outDir, basename)
            os.system(call)
