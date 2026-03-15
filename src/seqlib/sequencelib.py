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
    """Return the probability of a DNA sequence under a background GC model.

    Assumes each position is independently drawn from a 4-letter alphabet
    where G and C each have probability ``pGC/2`` and A and T each have
    probability ``(1-pGC)/2``.  The joint probability is the product of
    per-position probabilities.

    Args:
        seq: A DNA sequence string containing only ``A``, ``T``, ``G``,
            or ``C`` characters (upper case).
        pGC: The background GC probability in [0, 1].  Defaults to 0.5.

    Returns:
        The probability of observing ``seq`` under the model as a float.

    Raises:
        AssertionError: If ``pGC`` is outside [0, 1].
        ValueError: If ``seq`` contains a character other than
            ``A``, ``T``, ``G``, or ``C``.
    """
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
    """Transcribe a DNA sequence to RNA by replacing thymine with uracil.

    Performs a simple string substitution of every uppercase ``'T'``
    with ``'U'``.  Lower-case ``'t'`` characters are not converted.

    Args:
        seq: A DNA sequence string (upper case ``T`` will be replaced).

    Returns:
        The RNA sequence string with all ``'T'`` characters replaced by
        ``'U'``.
    """
    RNA = seq.replace('T', 'U')
    return RNA

def GenRandomSeq(length, type='DNA'):
    """Generate a random nucleotide sequence of a given length.

    Each position is drawn uniformly and independently from the
    four-letter alphabet appropriate for the requested sequence type.

    Args:
        length: The number of nucleotides in the returned sequence.
        type: The sequence type: ``'DNA'`` (alphabet ``A``, ``T``, ``G``,
            ``C``) or ``'RNA'`` (alphabet ``A``, ``U``, ``G``, ``C``).
            Defaults to ``'DNA'``.

    Returns:
        A random nucleotide sequence string of the specified length.
    """
    if type == 'DNA':
        chars = ['A','T','G','C']
    if type == 'RNA':
        chars = ['A','U','G','C']
    return ''.join([random.choice(chars) for i in range(length)])

def seed():
    """Re-seed the random number generator from the current system time.

    Calls :func:`random.seed` with no arguments, which uses the OS entropy
    source or the current time as the seed.  Useful for resetting
    deterministic state after a fixed seed has been set elsewhere.
    """
    random.seed()

def draw(distribution):
    """Draw a random index from a discrete probability distribution.

    Iterates through the distribution, accumulating a running sum, and
    returns the index of the first element whose cumulative sum exceeds
    a uniformly drawn random number.

    Args:
        distribution: A list of non-negative floats that sum to
            approximately 1.0.  Element ``i`` represents the probability
            of returning index ``i``.

    Returns:
        An integer index into ``distribution`` sampled according to the
        distribution's probabilities, or ``None`` if no element was
        selected (which can occur when probabilities do not sum to 1).
    """
    sum=0
    r = random.random()
    for i in range(0,len(distribution)):
        sum += distribution[i]
        if r< sum:
            return i

def makeDistFromFreqs(freqs):
    """Build a cumulative distribution list from a nucleotide frequency dict.

    Converts a dictionary of base frequencies into a list of cumulative
    boundary values suitable for use with :func:`draw`.  Bases are
    processed in the fixed order ``A``, ``T``, ``C``, ``G``.

    Args:
        freqs: A dictionary mapping nucleotide characters (``'A'``,
            ``'T'``, ``'C'``, ``'G'``) to their relative frequencies.
            Values should be non-negative and sum to 1.0.

    Returns:
        A list of five floats: the initial ``0.0`` followed by the
        cumulative sum after adding each of ``A``, ``T``, ``C``, ``G``
        in that order.
    """
    res = []
    chars = ['A','T','C','G']
    cum = 0
    res.append(cum)
    for i in chars:
        cum += freqs[i]
        res.append(cum)
    return res

def genRandomFromDist(length,freqs):
    """Generate a random DNA sequence drawn from a given base-frequency distribution.

    Builds a cumulative distribution from ``freqs`` and samples each
    position independently using :func:`draw`.

    Args:
        length: The number of nucleotides in the returned sequence.
        freqs: A dictionary mapping nucleotide characters (``'A'``,
            ``'T'``, ``'C'``, ``'G'``) to their probabilities.  Values
            should be non-negative and sum to 1.0.

    Returns:
        A random DNA sequence string of the specified length, with each
        base sampled proportionally to its frequency.
    """
    myDist = makeDistFromFreqs(freqs)
    chars = ['A','T','C','G']
    return ''.join([chars[draw(myDist)] for i in range(length)])

###########
#Motif Tools
###########
def allindices(string, sub, listindex=[], offset=0):
    """Find all start indices of substring ``sub`` within ``string``.

    Searches for non-overlapping occurrences of ``sub`` in ``string``
    starting from ``offset`` and appends each found index to
    ``listindex``.

    Warning:
        ``listindex`` uses a mutable default argument.  Repeated calls
        without explicitly passing a new list will accumulate results
        across calls.

    Args:
        string: The string to search within.
        sub: The substring to search for.
        listindex: A list to which found indices are appended.
            Defaults to a shared mutable list (see warning above).
        offset: The character position at which to start the search.
            Defaults to 0.

    Returns:
        The ``listindex`` list (same object passed in) with the start
        positions of all occurrences of ``sub`` appended.
    """
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex

def find_all(seq, sub):
    """Find all start positions of a substring within a sequence string.

    Iterates through ``seq`` looking for non-overlapping occurrences of
    ``sub`` using :func:`string.find` and collects each start index.

    Args:
        seq: The sequence string to search within.
        sub: The substring to search for.

    Returns:
        A list of integer start positions (0-based) of all occurrences of
        ``sub`` in ``seq``.  Returns an empty list if ``sub`` is not found.
    """
    #print "Looking for %s in %s"%(sub,seq)
    found = []
    next = string.find(seq,sub)
    while next != -1:
        found.append(next)
        next = string.find(seq,sub,next+1)
    return found

def kmer_dictionary_counts(seq,k,dic={}):
    """Count all k-mers in a sequence and store the counts in a dictionary.

    Slides a window of width ``k`` across ``seq`` and increments the
    count for each k-mer substring encountered.

    Warning:
        ``dic`` uses a mutable default argument.  Repeated calls without
        explicitly passing a fresh dict will accumulate counts across calls.

    Args:
        seq: The nucleotide (or any) sequence string to count k-mers in.
        k: The length of each k-mer.
        dic: A dictionary to update with k-mer counts.  Defaults to a
            shared mutable dict (see warning above).

    Returns:
        The updated ``dic`` dictionary mapping each k-mer string to its
        occurrence count in ``seq``.
    """
    for i in range(0, len(seq)-k):
        subseq = seq[i:][:k]
        #if not dic.has_key(subseq): dic[subseq] = 1
        #else: dic[subseq] = dic[subseq] + 1
        #OR
        dic[subseq] = 1 + dic.get(subseq,0)
    return dic

def kmer_dictionary(seq,k,dic={},offset=0):
    """Build a dictionary mapping each k-mer to its start positions in a sequence.

    Slides a window of width ``k`` across ``seq`` and records each
    1-based start position under the corresponding k-mer key.

    Warning:
        ``dic`` uses a mutable default argument.  Repeated calls without
        passing a fresh dict will accumulate positions across calls.

    Args:
        seq: The nucleotide (or any) sequence string to index.
        k: The length of each k-mer.
        dic: A dictionary to update with k-mer position lists.  Defaults
            to a shared mutable dict (see warning above).
        offset: Unused parameter retained for API compatibility.

    Returns:
        The updated ``dic`` dictionary mapping each k-mer string to a list
        of 1-based integer start positions at which it occurs in ``seq``.
    """
    for i in range(0,len(seq)-k):
        subseq = seq[i:][:k]
        dic.setdefault(subseq,[]).append(i+1)
    return dic

def kmer_stats(kmer,dic,genfreqs):
    """Compute enrichment statistics for a k-mer relative to genomic background.

    Calculates the signal-to-noise ratio (SNR) and Z-score for the
    observed count of ``kmer`` in a sequence compared to the count
    expected under a genomic-frequency background model.

    The expected count is ``sum(dic.values()) * genfreqs[kmer]``.

    Args:
        kmer: The k-mer string to evaluate.
        dic: A dictionary mapping k-mer strings to their observed counts,
            as returned by :func:`kmer_dictionary_counts`.
        genfreqs: A dictionary mapping k-mer strings to their expected
            background frequencies (floats summing to 1 across all k-mers
            of that length).

    Returns:
        A dict with keys ``'snr'`` (signal-to-noise ratio) and
        ``'zscore'`` (Z-score) if ``kmer`` is present in both ``dic`` and
        ``genfreqs``.  Returns ``None`` if ``dic`` is empty or ``kmer``
        is absent from either dictionary.
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
    """Collect and count 7-mer seeds from an iterable of sequence records.

    Iterates over sequence records, converts each from colorspace to DNA
    (by calling ``CSToDNA()`` on each record), extracts a 7-base seed
    from positions 1–7 (1-based, i.e. ``sequence[1:8]``), and counts the
    occurrences of each seed.  Prints progress every 10 000 records.

    Warning:
        ``seeds`` uses a mutable default argument.  Repeated calls
        without passing a fresh dict will accumulate counts across calls.

    Args:
        iter: An iterable of sequence-record objects.  Each object must
            have a ``sequence`` attribute and a ``CSToDNA()`` method that
            converts colorspace encoding to DNA in-place.
        seeds: A dictionary to update with seed counts.  Defaults to a
            shared mutable dict (see warning above).

    Returns:
        The updated ``seeds`` dictionary mapping 7-mer seed strings to
        their occurrence counts.
    """
    counter = 0
    for i in iter:
        counter+=1
        if counter%10000==0:
            print("%d" % counter)
        i.CSToDNA()
        seed = i.sequence[1:8]
        seeds[seed] = 1 + seeds.get(seed,0)
    return seeds
