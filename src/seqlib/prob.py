#!/usr/bin/env python
"""Probability and statistics tools for DNA sequence analysis.

Provides signal-to-noise ratio, Z-score, binning, cumulative sums,
nucleotide frequency utilities, Gaussian evaluation, moving averages,
Poisson and binomial probability functions, combinatorics, and
dictionary utility functions used throughout seqlib.
"""
import math
import operator
import random
import sys
from functools import reduce

import numpy as np


#######
#Probability Tools for DNA sequence analysis
#######
def snr(observed,expected):
    """Compute the signal-to-noise ratio (SNR) of an observed count vs an expected count.

    Calculates the simple ratio::

        SNR = observed / expected

    Args:
        observed: The observed count or value (numeric).
        expected: The expected count or value (numeric, must be non-zero).

    Returns:
        The ratio observed / expected as a float.
    """
    return observed/expected

def zscore(observed,expected):
    """Compute the Z-score of an observed count under a Poisson null model.

    Assumes the standard deviation equals the square root of the
    expected count (Poisson approximation)::

        Z = (observed - expected) / sqrt(expected)

    Args:
        observed: The observed count or value (numeric).
        expected: The expected count or value (numeric, must be positive).

    Returns:
        The Z-score as a float.
    """
    return (observed-expected)/math.sqrt(expected)

def which_bin(bins, x, safe=0):
    """Determine which bin interval a value ``x`` falls into.

    Given sorted bin boundary values, returns the 0-based index of the
    interval that contains ``x``.  For example, with boundaries
    ``[0, 5, 10, 15]``::

        x < 0       -> -1
        0 <= x < 5  ->  0
        5 <= x < 10 ->  1
        10 <= x < 15->  2
        x >= 15     ->  3   (or len(bins) when safe=0)

    Args:
        bins: A sorted list of numeric bin boundary values.
        x: The value to bin.
        safe: If ``1`` and ``x`` exactly equals ``bins[-1]``, returns
            ``len(bins)`` instead of the usual out-of-range value.
            Defaults to 0.

    Returns:
        An integer bin index.  Returns ``-1`` if ``x < bins[0]``, or
        ``len(bins)`` if ``x >= bins[-1]`` (unless ``safe=1`` applies).
    """
    if x<bins[0]: return -1
    for i in range(1,len(bins)):
        if x<bins[i]: return i-1
    if safe and x==bins[-1]: return len(bins)
    return len(bins)

def cumulative_sum(quality):
    """Compute the cumulative sum of a list in-place-style (returns a new list).

    Creates a copy of ``quality`` and then replaces each element with
    the running total up to and including that position.

    Args:
        quality: A list of numeric values.

    Returns:
        A new list of the same length as ``quality`` where element ``i``
        is the sum of ``quality[0]`` through ``quality[i]``.  Returns
        the input unchanged (empty list or falsy value) if ``quality``
        is empty.
    """
    if not quality: return quality
    sum_q = quality[:]
    for i in range(1,len(quality)):
        sum_q[i] = sum_q[i-1]+quality[i]
    return sum_q

def frequency_dic(seq):
    """Build a nucleotide frequency dictionary from a DNA sequence.

    Converts ``seq`` to uppercase and counts each of the four standard
    bases (A, C, G, T) as a fraction of the total sequence length.

    Args:
        seq: A DNA sequence string.  Mixed-case input is handled by
            uppercasing before counting.

    Returns:
        A dictionary mapping each of ``'A'``, ``'C'``, ``'G'``,
        ``'T'`` to its relative frequency (float in [0, 1]).  Bases not
        present in ``seq`` are mapped to 0.0.
    """
    dic = {}
    bases = ['A','C','G','T']
    seq=seq.upper()
    for b in bases:
        dic[b]=seq.count(b)/float(len(seq))
    return dic

def pick_one(dic):
    """Sample a single item from a dictionary of items and their probabilities.

    Builds a cumulative distribution from the dictionary's values and
    draws one item proportionally.  For example, with
    ``{'A': .18, 'C': .32, 'G': .32, 'T': .18}``, ``'A'`` is returned
    with probability 0.18, ``'C'`` with 0.32, and so on.

    Note:
        The function relies on :func:`cget` to extract values from
        ``dic.items()``; the behaviour may vary depending on dictionary
        iteration order (insertion order in Python 3.7+).

    Args:
        dic: A dictionary mapping hashable items to their relative
            probabilities.  Values should be non-negative; they need
            not sum to exactly 1.

    Returns:
        A randomly selected key from ``dic``, sampled with probability
        proportional to its value.
    """
    # {'A': .18, 'C': .32, 'G': .32, 'T': .18}
    # will generate A with probability .18 and so on
    items = dic.items()
    cums = cumulative_sum(cget(items,1))
    if 1: #debug:
        #print cums
        x = random.uniform(0,cums[-1])
        bin = which_bin(cums, x, safe=1)
        #print "%s is in bin %s and char %s. items=%s"%(
        #    x,bin,items[bin][0],items)
        return items[bin+1][0]
    else:
        return items[which_bin(cums, random.uniform(0,cums[-1]), safe=1)][0]

def pick_many(dic, n):
    """Sample ``n`` items independently from a probability dictionary.

    Builds a cumulative distribution from the dictionary's values once
    and then draws ``n`` samples with replacement.  For example, with
    ``{'A': .18, 'C': .32, 'G': .32, 'T': .18}``, each draw returns
    ``'A'`` with probability 0.18, ``'C'`` with 0.32, and so on.

    Args:
        dic: A dictionary mapping hashable items to their relative
            probabilities.  Values should be non-negative.
        n: The number of items to draw.

    Returns:
        A list of ``n`` keys from ``dic``, each sampled with probability
        proportional to its value.
    """
    # {'A': .18, 'C': .32, 'G': .32, 'T': .18}
    # will generate A with probability .18 and so on
    items = dic.items()
    cums = cumulative_sum(cget(items,1))
    choices = []
    for i in range(0,n):
        x = random.uniform(0,cums[-1])
        bin = which_bin(cums, x, safe=1)
        choices.append(items[bin+1][0])
    return choices

def gaussian(x,mu,sigma):
    """Evaluate the Gaussian (normal) PDF N(mu, sigma) at ``x``.

    Computes::

        f(x) = (1 / sqrt(2 * pi * sigma)) * exp(-((x - mu)^2) / (2 * sigma^2))

    Note:
        The normalisation constant uses ``sqrt(2*pi*sigma)`` rather than
        the more common ``sigma*sqrt(2*pi)``.  For the function to
        integrate to 1 the usual convention is ``sigma`` (standard
        deviation) in the denominator as ``sigma * sqrt(2*pi)``.

    Args:
        x: The point at which to evaluate the PDF.
        mu: The mean of the Gaussian.
        sigma: The standard deviation of the Gaussian.

    Returns:
        The PDF value at ``x`` as a float.
    """
    return ( (1.0/math.sqrt(2*math.pi*sigma)) * (math.e**(-((x-mu)**2)/(2*sigma**2))))

def make_gaussian(mu,sigma):
    """Create a Gaussian PDF function with fixed mean and standard deviation.

    Returns a callable that evaluates the Gaussian PDF at any point
    ``x``, with ``mu`` and ``sigma`` captured by closure.

    Example::

        N2_3 = make_gaussian(2, 3)
        N2_3(4)  # -> gaussian(4, mu=2, sigma=3)

    Args:
        mu: The mean of the Gaussian.
        sigma: The standard deviation of the Gaussian.

    Returns:
        A function ``f(x)`` that evaluates the Gaussian N(mu, sigma) at
        ``x``.
    """
    return lambda x,mu=mu,sigma=sigma: ( (1.0/math.sqrt(2*math.pi*sigma)) * (math.e**(-((x-mu)**2)/(2*sigma**2))))

def make_adder(n):
    """Create an adder function that adds a fixed value ``n`` to its argument.

    Returns a callable that adds ``n`` (captured by closure) to any
    input ``x``.

    Example::

        Add2 = make_adder(2)
        Add2(3)  # -> 5

    Args:
        n: The fixed value to add.

    Returns:
        A function ``f(x)`` that returns ``x + n``.
    """
    return lambda x,n=n: x+n

#############
#Math Primitives
#############
loge_2 = math.log(2)

def avg(l,precise=0):
    """Compute the arithmetic mean of a list of numbers.

    Args:
        l: A list of numeric values.
        precise: If non-zero, divide by ``float(len(l))`` for a
            floating-point result.  If 0 (default), divide by
            ``len(l)`` using integer or floor division.

    Returns:
        The mean of ``l`` as a number, or 0 if ``l`` is empty.
    """
    if not l: return 0
    if precise:
        return reduce(operator.add,l,0)/float(len(l))
    else:
        return reduce(operator.add,l,0)/len(l)

def movavg(s, n):
    """Compute an n-period moving average for a time series.

    Uses cumulative sums for an O(len(s)) implementation::

        MA[i] = mean(s[i-n+1 : i+1])

    The result has length ``len(s) - n + 1``.

    Args:
        s: A list or array of numeric values ordered from oldest
            (index 0) to most recent (index -1).
        n: The window size (number of periods) for the moving average.

    Returns:
        A NumPy array of the moving average values.  The array has
        ``len(s) - n + 1`` elements.
    """
    s = np.array(s)
    c = np.cumsum(s)
    return (c[n-1:] - c[:-n+1]) / float(n)


def median(l):
    """Compute the median of a list of numbers.

    Sorts ``l`` and returns the middle value for odd-length lists or the
    average of the two middle values for even-length lists.

    Args:
        l: A list of numeric values.

    Returns:
        The median value, or ``None`` if ``l`` is empty.
    """
    if not l: return None
    l = sorted(l)
    if len(l)%2: return sorted(l)[len(l)//2]
    else: return (l[len(l)//2]+l[len(l)//2-1])/2.0

def stdev(l, failfast=1):
    return math.sqrt(variance(l,failfast=failfast))

def variance(l,failfast=1):
    if (not l) or len(l)==1:
        if failfast: raise "tools.variance: Not enough samples.  Need >= 2, got %s"%len(l)
        else: return 0#'N/A'
    m = avg(l,1)
    s = 0
    for i in l:
        s = s + (i-m)*(i-m)
    return s / (len(l)-1)

def log2(x):
    #converting bases: log_a(b) = log_c(b)/log_c(a)
    #i.e. log_2(x) = log_e(2)/log_e(x) = log_10(2)/log_10(x)
    return math.log(x)/float(loge_2)

def log_k(x,k):
    return math.log(x)/math.log(k)

def prob2score(prob):
    #1/100 -> 20
    try:
        return -10*float(math.log10(float(prob)))
    except:
        return -1

def p2bits(p):
    """Takes p-value and returns negative log2"""
    return -log2(p)

def factorial(n):
    result = 1
    for i in range(n,0,-1):
        #print i
        result = result * i
    return result

###########
#Poisson
###########
def poisson_expected(rate):
    for x in range(1,50,1):
        p = poisson(rate,x)
        print(f"{x}\t{p}\t{12000000*p}")

def poisson(rate, x):
    """Returns the probability of observing a count of x"""
    return math.exp(-rate)*(rate**x)/factorial(x)

######################
#Binomial Distribution
#######################
def binomial_likelihood_ratio(ps,k,n):
    # p[0] is the null hypothesis
    # p[1] is the hypothesis being tested
    assert(len(ps)==2)
    likelihoods = []
    for p in ps:
        likelihoods.append(binomial(p,k,n))
    #i = argmax(likelihoods)
    #p = likelihoods[i] / sum(likelihoods)
    #return p
    if likelihoods[0]: return np.log(likelihoods[1]) / likelihoods[0]
    else:
        print("Warning: likelihood ratio set to sys.maxsize.  p(H1)=%s, p(H0)=0"%(p[1]))
        return sys.maxsize

def binomial_log_likelihood_ratio(ps,k,n):
    return log_binomial(ps[1],k,n) - log_binomial(ps[0],k,n)

def log_binomial(p,k,n):
    # the log probability of seeing exactly k successes in n trials
    # given the probability of success is p
    return log_n_choose_k(n,k)+math.log(p)*k+math.log(1-p)*(n-k)

def binomial(p,k,n):
    # probability of seeing exactly k successes in n trials, given
    # the probability of success is p
    #return n_choose_k(n,k)*(p**k)*((1-p)**(n-k))
    return n_choose_k(n,k)*(p**k)*((1-p)**(n-k))

def cumBinomial(p,k,n):
    #Returns the cumulative probability from the binomaial distribution
    Pval = 0.0
    for j in range(0,k+1):
        Pval+=binomial(p,j,n)
    return Pval

def n_choose_k(n,k):
    # (n k) = n! / (k! (n-k)!)
    #
    #         n*(n-1)*(n-2)*....*(n-k+1)
    #       = --------------------------
    #              k*(k-1)*...*1
    assert(k<=n)
    k = min(k, n-k)
    nominator   = range(n,n-k,-1)
    denominator = range(k,0,-1)

    result = 1.0
    for nom, den in zip(nominator, denominator):
        result = (result * nom) / den
        #result = result*nom
        #print result
        #result = result/den
        #print result

    return result

def log_n_choose_k(n,k):
    # (n k) = n! / (k! (n-k)!)
    #
    #         n*(n-1)*(n-2)*....*(n-k+1)
    #       = --------------------------
    #              k*(k-1)*...*1
    assert(k<=n)
    k = min(k, n-k)
    nominator   = range(n,n-k,-1)
    denominator = range(k,0,-1)

    result = 0
    for nom, den in zip(nominator, denominator):
        result = (result + math.log(nom)) - math.log(den)
    return result

#################
#Dictionary Tools
#################
def cget(diclist, key, strict=1):
    # cross_get was: gather(diclist,key)
    # gathers the same key from a list of dictionaries
    # can also be used in lists

    # input: a list of dictionaries all of which contains key
    # output: a list of elements d[key] for each d in diclist
    if strict:
        # return map(lambda d,key=key: d[key], diclist)
        result = [None]*len(diclist)
        for i in range(0,len(diclist)):
            result[i] = diclist[i][key]
        return result
    else:
        results = []
        for dic in diclist:
            if dic and generic_has_key(dic,key):
                results.append(dic[key])
        return results
