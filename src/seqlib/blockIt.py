'''
Block-iT miRNA expression vector insert design utilities.

Given a 21-mer siRNA candidate sequence, generates the forward and reverse
oligonucleotide sequences required for cloning into the pcDNA6.2-GW/Em-GFP/miR
expression vector (Invitrogen Block-iT Kit).

Created on Oct 14, 2009
Takes as input a 21mer sequence (candidate siRNA) and creates the appropriate fwd and rev oligo
sequences to order for insertion into the pcDNA6.2-GW/Em-GFP/miR expression vector from
Invitrogen (Block-It Kit)

@author: lgoff
'''
import sys

from . import sequencelib as sequence

fwdAdapter = 'TGCTG'
loopSequence = 'GTTTTGGCCACTGACTGAC'
revAdapter = 'CCTG'

def makeBlockItInsert(seq):
    """Design forward and reverse oligos for a Block-iT miRNA insert.

    Constructs the forward strand by concatenating the fixed forward adapter,
    the reverse complement of seq, the loop sequence, and a modified copy of
    seq (positions 0-7 joined directly to positions 10 onward, skipping 8-9).
    The reverse strand is the reverse complement of the forward strand
    (excluding the first four adapter bases), prefixed by the reverse adapter.

    Args:
        seq: A 21-nucleotide DNA string representing the candidate siRNA
            sense sequence (5' to 3').

    Returns:
        A tuple (fwdStrand, revStrand) where both elements are DNA strings
        suitable for ordering as oligonucleotides.
    """
    fwdStrand = fwdAdapter+sequence.reverse_complement(seq)+loopSequence+seq[:8]+seq[10:]
    revStrand = revAdapter+sequence.reverse_complement(fwdStrand[4:])
    return (fwdStrand, revStrand)

def printBlockIt(seqs):
    """Print the forward and reverse oligo sequences from a Block-iT insert tuple.

    Prints each strand labeled 'FWD' or 'REV' to stdout. Also computes a
    base-pairing alignment string between the forward and reverse strands,
    though the alignment is computed but not printed.

    Args:
        seqs: A tuple (fwdStrand, revStrand) as returned by makeBlockItInsert.
    """
    print("FWD:\t%s" % seqs[0])
    print("REV:\t%s" % seqs[1])

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
