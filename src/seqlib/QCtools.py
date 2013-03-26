#!/usr/bin/env python
'''
Created on May 6, 2010

@author: lgoff
'''
import numpy as np
import re



def makePWM(fastqFile,readLen,freq=True):
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
        for j in xrange(0,len(i['sequence'])):
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
        if line[0] <> "@":
            raise ValueError("Records in csfastq files should start with '@'")
        name = line[1:].rstrip()
        line = handle.readline()
        sequence = line.rstrip()
        line = handle.readline()
        if line[0] <> "+":
            raise ValueError("Fastq file does not appear to be formatted correctly")
        line = handle.readline()
        quals = line.rstrip()
        fastq = {'name':name,'sequence':sequence,'quals':quals}
        yield fastq
        line = handle.readline()
        if not line: return
    assert False, "Should not reach this line"
