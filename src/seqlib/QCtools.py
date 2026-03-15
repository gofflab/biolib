#!/usr/bin/env python
'''
Quality control tools for sequencing data.

Provides a FASTQ file parser and a position-weight matrix (PWM) builder for
inspecting base-composition biases across read positions.

Created on May 6, 2010

@author: lgoff
'''

import numpy as np


def makePWM(fastqFile,readLen,freq=True):
    """Build a position-weight matrix of base composition from a FASTQ file.

    Iterates over all records in a FASTQ file and tallies the occurrence of
    each nucleotide (A, C, G, T) at every position across ``readLen``
    positions.  Ambiguous bases (e.g. 'N') are silently ignored.
    Optionally converts raw counts to per-position frequencies.

    Args:
        fastqFile: Path to the FASTQ file to process.
        readLen: Expected read length (number of positions to track).
        freq: If True (default), each base count vector is divided by the
            total count at that position to produce a frequency.  If False,
            raw counts are returned.

    Returns:
        A dict with keys 'A', 'C', 'G', 'T', and 'Total'.  Each key maps to
        a numpy array of length ``readLen``.  The 'Total' array contains the
        total number of valid base observations at each position; the
        individual base arrays contain either counts or frequencies depending
        on the ``freq`` argument.
    """
    bases = ['A','C','G','T']
    pwm = {
           'A':np.zeros(readLen),
           'C':np.zeros(readLen),
           'G':np.zeros(readLen),
           'T':np.zeros(readLen),
           'Total':np.zeros(readLen)
           }


    #Iterate over fastq records
    iter=FastqIterator(fastqFile)
    for i in iter:
        for j in range(0,len(i['sequence'])):
            try:
                pwm[i['sequence'][j]][j] += 1
                pwm['Total'][j] += 1
            except KeyError:
                pass
    if freq:
        for base in bases:
            pwm[base] = pwm[base]/pwm['Total']
    return pwm

################
#Parsers
################
def FastqIterator(fastqFile):
    """Iterate over records in a FASTQ file.

    Skips any non-FASTQ header text at the start of the file (lines that do
    not begin with '@') and then yields one dict per record.  The file is
    expected to use standard four-line FASTQ format: a '@'-prefixed name
    line, a sequence line, a '+' line, and a quality line.

    Args:
        fastqFile: Path to the FASTQ file to parse.

    Yields:
        A dict with keys:
            ``'name'``: Read name string (the '@' prefix is stripped).
            ``'sequence'``: Nucleotide sequence string.
            ``'quals'``: ASCII quality string.

    Raises:
        ValueError: If a record's name line does not start with '@'.
        ValueError: If the separator line between sequence and qualities
            does not start with '+'.
    """
    handle = open(fastqFile,'r')
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "": return
        if line [0] == "@":
            break

    #Begin walk through csfasta records
    while True:
        if not line: break
        if line[0] != "@":
            raise ValueError("Records in csfastq files should start with '@'")
        name = line[1:].rstrip()
        line = handle.readline()
        sequence = line.rstrip()
        line = handle.readline()
        if line[0] != "+":
            raise ValueError("Fastq file does not appear to be formatted correctly")
        line = handle.readline()
        quals = line.rstrip()
        fastq = {'name':name,'sequence':sequence,'quals':quals}
        yield fastq
        line = handle.readline()
        if not line: return
    assert False, "Should not reach this line"
