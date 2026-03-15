#/usr/bin/env python
"""Sequence utility functions for DNA/RNA analysis.

Provides parsers, generic sequence tools, and motif tools for working
with biological sequence data including FASTA parsing, complement
computation, GC content, k-mer analysis, and random sequence generation.
"""
import math
import operator
import random
import string

from . import prob


######
#Parsers
######
def FastaIterator(handle):
    """Iterate over FASTA records in an open file handle.

    Skips any header text before the first '>' character, then yields
    one record dict per FASTA entry.  Each sequence has internal
    whitespace stripped and lines joined into a single string.

    Args:
        handle: A readable file object (e.g. opened with ``open(path, 'r')``)
            positioned at or before the first FASTA record.

    Yields:
        A dict with keys:
            ``'name'``: The record header string (everything after ``>``
            on the header line, whitespace-stripped).
            ``'sequence'``: The concatenated sequence string with all
            internal spaces removed.

    Raises:
        ValueError: If a record block does not begin with a ``>``
            character as required by the FASTA format.
    """
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break

    while True:
        if line[0] !=">":
            raise ValueError("Records in Fasta files should start with a '>' character")
        name = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            lines.append(line.rstrip().replace(" ",""))
            line = handle.readline()
        #Return record then continue
        newSeq = {'name':name,'sequence':"".join(lines)}
        yield newSeq

        if not line : return #StopIteration
    assert False, "Should not reach this line"

bed_fields = ['chr','start','end','label','score','strand']

###
#Generic Sequence tools
###

def complement(s):
    """Return the base-by-base complement of a DNA sequence as a list.

    Handles both upper- and lower-case input characters.  Note that the
    lower-case mapping contains a known quirk: ``'c'`` maps to ``'t'``
    instead of ``'g'``.

    Args:
        s: An iterable of single-character DNA bases (``A``, ``T``, ``G``,
            ``C`` in either case).

    Returns:
        A list of single-character strings representing the complemented
        bases in the same order as the input.

    Raises:
        KeyError: If a character in ``s`` is not present in the complement
            lookup table.
    """
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'a': 't', 'c': 't', 'g': 'c', 't': 'a'
            }
    complseq = [comp[base] for base in s]
    return complseq

def reverse_complement(s):
    """Return the reverse complement of a DNA sequence string.

    Reverses the sequence and then complements each base using
    :func:`complement`.

    Args:
        s: A DNA sequence string containing bases ``A``, ``T``, ``G``, ``C``
            (upper or lower case).

    Returns:
        A string that is the reverse complement of ``s``.
    """
    seq = list(s)
    seq.reverse()
    return ''.join(complement(seq))

def rcomp(s):
    """Return the reverse complement of an uppercase DNA string.

    Uses ``str.translate`` with a precomputed translation table for
    ``A<->T`` and ``C<->G``, then reverses the result with a slice.
    Equivalent to :func:`reverse_complement` but operates only on
    uppercase bases via the translation table.

    Args:
        s: An uppercase DNA sequence string (``A``, ``T``, ``C``, ``G``).

    Returns:
        A string that is the reverse complement of ``s``.
    """
    return s.translate(string.maketrans("ATCG","TAGC"))[::-1]

def getTm(seq):
    """Calculate the melting temperature (Tm) of a DNA sequence.

    Uses the nearest-neighbour-inspired empirical formula::

        Tm = 79.8 + 18.5*log10([Na+]) + 58.4*GC + 11.8*GC^2 - 820/len

    where ``[Na+]`` is fixed at 0.05 M and ``GC`` is the fractional GC
    content of the sequence.

    Args:
        seq: A DNA sequence string.

    Returns:
        The estimated melting temperature in degrees Celsius as a float.
    """
    Tm = 79.8 + 18.5*math.log10(0.05) + (58.4 * getGC(seq)) + (11.8 * getGC(seq)**2) - (820/len(seq))
    return Tm

def getGC(seq):
    """Return the fractional GC content of a DNA sequence.

    Counts both upper- and lower-case ``G`` and ``C`` characters and
    divides by the total sequence length.

    Args:
        seq: A DNA sequence string.

    Returns:
        A float in [0, 1] representing the proportion of G and C bases.
    """
    return (seq.count('C')+seq.count('G')+seq.count('c')+seq.count('g'))/float(len(seq))

def gc_content(seq):
    """Return the percentage GC content of a nucleotide sequence.

    Counts G and C characters (upper and lower case) and divides by the
    sum of all A, T, U, G, C characters (upper and lower case), ignoring
    any ambiguity codes or gap characters.  The result is scaled to
    percentage (0–100).

    Args:
        seq: A DNA or RNA sequence string.

    Returns:
        A float representing GC content as a percentage (0–100).
    """
    gc = mcount(seq, 'GCgc')
    at = mcount(seq, 'ATUatu')
    return 100*gc/float((gc+at))

def mcount(s, chars):
    """Count all occurrences of any character in ``chars`` within string ``s``.

    Iterates over each character in ``chars`` and accumulates the count of
    its appearances in ``s`` using ``string.count``.

    Args:
        s: The string to search within.
        chars: A string whose individual characters are each counted in ``s``.

    Returns:
        The total number of occurrences of any character from ``chars``
        found in ``s`` as an integer.
    """
    # sums the counts of appearances of each char in chars
    count = 0
    for char in chars:
        count = count+string.count(s,char)
    return count

def prob_seq(seq, pGC=.5):
    # given a GC content, what is the probability
    # of getting the particular sequence

    assert(0<=pGC<=1)
    # the probability of obtaining sequence seq
    # given a background gc probability of .5
    ps = []
    for char in seq:
        if char in 'CG': ps.append(pGC/2)
        elif char in 'AT': ps.append((1-pGC)/2)
        else: raise ValueError("Unexpected char: " + repr(char))
    return reduce(operator.mul, ps, 1)

def transcribe(seq):
    RNA = seq.replace('T', 'U')
    return RNA

def GenRandomSeq(length, type='DNA'):
    if type == 'DNA':
        chars = ['A','T','G','C']
    if type == 'RNA':
        chars = ['A','U','G','C']
    return ''.join([random.choice(chars) for i in range(length)])

def seed():
    random.seed()

def draw(distribution):
    sum=0
    r = random.random()
    for i in range(0,len(distribution)):
        sum += distribution[i]
        if r< sum:
            return i

def makeDistFromFreqs(freqs):
    res = []
    chars = ['A','T','C','G']
    cum = 0
    res.append(cum)
    for i in chars:
        cum += freqs[i]
        res.append(cum)
    return res

def genRandomFromDist(length,freqs):
    """Generates a random sequence of length 'length' drawing from a distribution of
    base frequencies in a dictionary"""
    myDist = makeDistFromFreqs(freqs)
    chars = ['A','T','C','G']
    return ''.join([chars[draw(myDist)] for i in range(length)])

###########
#Motif Tools
###########
def allindices(string, sub, listindex=[], offset=0):
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex

def find_all(seq, sub):
    #print "Looking for %s in %s"%(sub,seq)
    found = []
    next = string.find(seq,sub)
    while next != -1:
        found.append(next)
        next = string.find(seq,sub,next+1)
    return found

def kmer_dictionary_counts(seq,k,dic={}):
    """Returns a dictionary of k,v = kmer:'count of kmer in seq'"""
    for i in range(0, len(seq)-k):
        subseq = seq[i:][:k]
        #if not dic.has_key(subseq): dic[subseq] = 1
        #else: dic[subseq] = dic[subseq] + 1
        #OR
        dic[subseq] = 1 + dic.get(subseq,0)
    return dic

def kmer_dictionary(seq,k,dic={},offset=0):
    """Returns dictionary of k,v = kmer:'list of kmer start positions in seq' """
    for i in range(0,len(seq)-k):
        subseq = seq[i:][:k]
        dic.setdefault(subseq,[]).append(i+1)
    return dic

def kmer_stats(kmer,dic,genfreqs):
    """Takes as argument a kmer string, a dictionary with kmers as keys from kmer_dictionary_counts, and a dictionary
        of genomic frequencies with kmers as keys. Returns a dictionary of stats for kmer ("Signal2Noise Ratio, Z-score")
    """
    if not dic: return
    if kmer in dic.keys() and kmer in genfreqs.keys():
        observed = dic[kmer]
        expected = sum(dic.values())*genfreqs[kmer]
        snr = prob.snr(observed,expected)
        zscore = prob.zscore(observed, expected)
        return {'snr':snr,'zscore':zscore}
    else: return

def get_seeds(iter,seeds={}):
    counter = 0
    for i in iter:
        counter+=1
        if counter%10000==0:
            print("%d" % counter)
        i.CSToDNA()
        seed = i.sequence[1:8]
        seeds[seed] = 1 + seeds.get(seed,0)
    return seeds
