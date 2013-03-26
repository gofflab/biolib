#!/usr/bin/env python
'''
Created on Oct 8, 2009
Generates list of candidate siRNAs from .fasta sequence given as argument

@author: lgoff
'''

"""
http://www.protocol-online.org/prot/Protocols/Rules-of-siRNA-design-for-RNA-interference--RNAi--3210.html
"""
import sequencelib
import math,sys,blockIt
    
def main(fastaFile):
    """Do it all"""
    handle = open(fastaFile,'r')
    iter = sequencelib.FastaIterator(handle)
    for i in iter:
        print "%s|Candidate siRNAs:" % (i['name'])
        evaluateSequence(i["sequence"])
        
def evaluateSequence(seq,scoreCutoff=6):
    """Wrapper for testCandidate() that iterates across sequence provided and returns candidates with a score >= scoreCutoff (default = 6)"""
    for i in range(0,len(seq)-21):
        candidate = seq[i:i+21]
        score = testCandidate(candidate)
        if score>=6:
            print "\t%d\t%s\t%.2f" % (i,candidate,score),
            insertSeqs = blockIt.makeBlockItInsert(candidate)
            print "Fwd:%s\tRev:%s" % (insertSeqs[0],insertSeqs[1]) 
            
def testCandidate(seq):
    """Checks 21mer candidates against siRNA rules and assigns a score on a scale of 0-8"""
    #seq = seq.upper()
    if len(seq)!=21:
        assert ValueError("Candidate is not 21nt in length")
        return False
    score = 0.0
    gc = getGC(seq)
    #Criteria 1: Moderate to low (30%-52%) GC Content (1 point)
    if 0.3 >= gc and gc <= 0.52:
        score += 1
    #Criteria 2: At least 3 A/Us at positions 15-19 (sense) (1 point /per A or U)
    tmp = seq[14:18].count('A')+seq[14:18].count('T')+seq[14:18].count('t')+seq[14:18].count('a')
    if tmp>=3:
        score += tmp
    #Criteria 3: Lack of internal repeats (Tm<20 degrees C) (1 point)
    Tm = getTm(seq)
    if Tm<20.0:
        score += 1
    #Criteria 4: A at position 19 (sense) (1 point)
    if seq[18] in ['A','a']:
        score += 1
    #Criteria 5: A at position 3 (sense) (1 point)
    if seq[2] in ['A','a']:
        score += 1
    #Criteria 6: U at position 10 (sense) (1 point)
    if seq[9] in ['T','t']:
        score += 1
    #Criteria 7: No G/C at position 19 (sense) (-1 point)
    if seq[18] in ['G','g'] or seq[18] in ['C','c']:
        score -= 1
    #Criteria 8: No G at position 13 (sense) (-1 point)
    if seq[12] in ['G','g']:
        score -= 1
    #Criteria 9: No stretches of 4 or more bases (-5 point)
    for i in ['A','C','G','T','a','c','g','t']:
        if seq.count(i*4)>0:
            score -= 5
    return score

def getTm(seq):
    Tm = 79.8 + 18.5*math.log10(0.05) + (58.4 * getGC(seq)) + (11.8 * getGC(seq)**2) - (820/len(seq))
    return Tm

def getGC(seq):
    seq = seq.upper()
    return (seq.count('C')+seq.count('G'))/float(len(seq))

######
#dsRNA rules from Vera et al. (updated 2-1-10)
######
def scanPromoter(promSeq):
    """
    Evaluates candidate dsRNAs for RNAa from a given sequence.  Returns a list of dictionaries of candidates and their score.
    """
    promSeq = promSeq.upper()
    window = 19
    candidates = []
    
    for i in range(len(promSeq)-window):
        candidates.append({})
        candidates[i]['seq'] = promSeq[i:i+window]
        candidates[i]['pos'] = -(len(promSeq)-i)
        candidates[i]['gc'] = getGC(candidates[i]['seq'])
        candidates[i]['score'] = 0.0
        
        #dsRNA Design Rules
        
        #GC content must be between 40-65%
        if 0.4 <= candidates[i]['gc'] and candidates[i]['gc'] <=0.65:
            candidates[i]['score'] += 1
        
        #Consecutive nucleotides >=4 are penalized
        for n in ['A','C','G','T','a','c','g','t']:
            if candidates[i]['seq'].count(n*4)>0:
                candidates[i]['score'] -= 5
            
        #19th position should be an 'A'
        if candidates[i]['seq'][18] in ['A','a']:
            candidates[i]['score'] += 1
            
        #Criteria 7: No G/C at position 19 (sense) (-1 point)
        if candidates[i]['seq'][18] in ['G','g'] or candidates[i]['seq'][18] in ['C','c']:
            candidates[i]['score'] -= 1
        
        #Position 18 should be an 'A' or 'T' preferrably an 'A'
        if candidates[i]['seq'][17] in ['A','a','T','t']:
            if candidates[i]['seq'][17] in ['A','a']:
                candidates[i]['score'] += 2
            if candidates[i]['seq'][17] in ['T','t']:
                candidates[i]['score'] += 1
        
        #Position 7 should be a 'T'
        if candidates[i]['seq'] in ['T','t']:
            candidates[i]['score'] += 1
        
        #The 20th-23rd positions (flanking the 3' end of a target) were preferably 'A's or 'T's
        tmp = promSeq[i+20:i+23].count('A')+promSeq[i+20:i+23].count('T')+promSeq[i+20:i+23].count('a')+promSeq[i+20:i+23].count('t')
        if tmp>=3:
            candidates[i]['score'] += tmp
        
        #Score for lack of internal repeats
        candidates[i]['Tm'] = getTm(candidates[i]['seq'])
        if candidates[i]['Tm']<20.0:
            candidates[i]['score'] += 1
            
    #Sort list by score
    return sorted(candidates,key=lambda k: k['score'],reverse=True)

def ASOscan(targetSeq):
    """
    Evaluates candidate dsRNAs for RNAa from a given sequence.  Returns a list of dictionaries of candidates and their score.
    """
    targetSeq = sequencelib.rcomp(targetSeq)
    window = 20
    candidates = []
    
    for i in range(len(targetSeq)-window):
        candidates.append({})
        candidates[i]['seq'] = targetSeq[i:i+window]
        candidates[i]['pos'] = -(len(targetSeq)-i)
        candidates[i]['gc'] = getGC(candidates[i]['seq'])
        candidates[i]['score'] = 0.0
        
        #dsRNA Design Rules
        
        #GC content must be between 40-65%
        if 0.45 <= candidates[i]['gc'] and candidates[i]['gc'] <=0.65:
            candidates[i]['score'] += 2
        
        #Consecutive nucleotides >=4 are penalized
        for n in ['A','C','G','T','a','c','g','t']:
            if candidates[i]['seq'].count(n*4)>0:
                candidates[i]['score'] -= 5
            
        #19th position should be an 'A'
        if candidates[i]['seq'][18] in ['A','a']:
            candidates[i]['score'] += 0
            
        #Criteria 7: No G/C at position 19 (sense) (-1 point)
        if candidates[i]['seq'][18] in ['G','g'] or candidates[i]['seq'][18] in ['C','c']:
            candidates[i]['score'] -= 0
        
        #Position 18 should be an 'A' or 'T' preferrably an 'A'
        if candidates[i]['seq'][17] in ['A','a','T','t']:
            if candidates[i]['seq'][17] in ['A','a']:
                candidates[i]['score'] += 0
            if candidates[i]['seq'][17] in ['T','t']:
                candidates[i]['score'] += 0
        
        #Position 7 should be a 'T'
        if candidates[i]['seq'] in ['T','t']:
            candidates[i]['score'] += 0
        
        #The 20th-23rd positions (flanking the 3' end of a target) were preferably 'A's or 'T's
        tmp = targetSeq[i+20:i+23].count('A')+targetSeq[i+20:i+23].count('T')+targetSeq[i+20:i+23].count('a')+targetSeq[i+20:i+23].count('t')
        if tmp>=3:
            #candidates[i]['score'] += tmp
            candidates[i]['score'] += 0
            
        #Score for lack of internal repeats
        candidates[i]['Tm'] = getTm(candidates[i]['seq'])
        if candidates[i]['Tm']>45.0:
            candidates[i]['score'] += 2
            
    #Sort list by score
    return sorted(candidates,key=lambda k: k['score'],reverse=True)

def makeDsRNA(seq):
    if len(seq)!=19:
        assert ValueError("Candidate is not 19nt in length")
        return False
    seq = seq.upper()
    revSeq = sequencelib.rcomp(seq)
    return ["r"+"r".join(seq)+"TT","r"+"r".join(revSeq)+"TT"]
                                         
def veraMain(fastaFile):
    """Do it all"""
    handle = open(fastaFile,'r')
    iter = sequencelib.FastaIterator(handle)
    for i in iter:
        print "-----------------------------------------------------------------\n%s Promoter Candidate dsRNAs\n-----------------------------------------------------------------" % (i['name'])
        candidates = scanPromoter(i['sequence'])
        for candidate in candidates[:10]:
            dsRNA = makeDsRNA(candidate['seq'])
            print "Pos:\t%d\nCandidate:\t%s\nScore:\t%.2f\nTm:\t%.2f\nGC:\t%.2f\nFwd:\t%s\nRev:\t%s\n------------------------" % (candidate['pos'],candidate['seq'],candidate['score'],candidate['Tm'],candidate['gc'],dsRNA[0],dsRNA[1])

def ASOMain(fastafile):
    """Takes a fasta sequnce of RNAs, reverse-complements and scans for ASO sequences"""
    handle = open(fastafile,'r')
    iter = sequencelib.FastaIterator(handle)
    for i in iter:
        print "----------------------------------------------------------\n%s ASO Candidate Regions (sequence is transcript-strand)\n---------------------------------------------------------" % (i['name'])
        candidates = ASOscan(i['sequence'])
        for candidate in candidates[:10]:
            #dsRNA = makeDsRNA(candidate['seq'])
            if candidate['seq'].count('a')+candidate['seq'].count('t')+candidate['seq'].count('g')+candidate['seq'].count('c') >0:
                continue
            else:
                print "Pos:\t%d\nCandidate:\t%s\nScore:\t%.2f\nTm:\t%.2f\nGC:\t%.2f\n------------------------" % (candidate['pos'],candidate['seq'],candidate['score'],candidate['Tm'],candidate['gc'])


if __name__=="__main__":
    VeraMain(sys.argv[1])