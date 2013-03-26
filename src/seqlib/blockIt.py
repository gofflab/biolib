'''
Created on Oct 14, 2009
Takes as input a 21mer sequence (candidate siRNA) and creates the appropriate fwd and rev oligo
sequences to order for insertion into the pcDNA6.2-GW/Em-GFP/miR expression vector from
Invitrogen (Block-It Kit)

@author: lgoff
'''
import sys
import sequencelib as sequence

fwdAdapter = 'TGCTG'
loopSequence = 'GTTTTGGCCACTGACTGAC'
revAdapter = 'CCTG'

def makeBlockItInsert(seq):
    fwdStrand = fwdAdapter+sequence.reverse_complement(seq)+loopSequence+seq[:8]+seq[10:]
    revStrand = revAdapter+sequence.reverse_complement(fwdStrand[4:])
    return (fwdStrand,revStrand)

def printBlockIt(seqs):
    """Takes as input the tuple returned from makeBlockItInsert and prints the result to stdout"""
    print "FWD:\t%s" % seqs[0]
    print "REV:\t%s" % seqs[1]
    
    alignment = '    '
    revRev = seqs[1][::-1]
    for i in range(len(seqs[1])-4):
        #print "%s%s" % (fwdStrand[i+4],sequence.reverse_complement(revRev[i]))
        if seqs[0][i+4]==sequence.reverse_complement(revRev[i]):
            alignment+="|"
        else:
            alignment+=" "
###
#Main
###    
if __name__ == '__main__':
    seq = sys.argv[1]
    makeBlockItInsert(seq)
    pass