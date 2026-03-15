"""Statistical and mathematical utilities for biological data analysis.

Provides descriptive statistics, probability distributions (PDF and CDF),
random variates, regression, sliding-window operations, curve fitting, and
special mathematical functions.  Functions that require external tools (rpy2,
gnuplot) fall back gracefully or raise ``NotImplementedError`` when those
dependencies are absent.
"""
# python libs
import cmath
import os
import random
from collections import Counter, defaultdict
from math import *

import numpy as np
import pandas as pd

# rasmus libs replaced with local imports and inlined utilities
# from rasmus import util       # removed: rasmus not Python 3 compatible
# from rasmus import algorithms # removed: use local algorithms module
# from rasmus import tablelib   # removed: replaced with pandas DataFrame
from . import algorithms


def prod(lst):
    """Compute the product of a list of positive numbers via log-space summation.

    Calculates ``exp(sum(log(i) for i in lst))``, which avoids numerical
    overflow for large lists by working in log space.  All values in
    ``lst`` must be strictly positive.

    Args:
        lst: An iterable of strictly positive numbers.

    Returns:
        The product of all elements in ``lst`` as a float.
    """
    return exp(sum(log(i) for i in lst))

def mean(vals):
    """Compute the arithmetic mean of a sequence of numbers.

    Iterates through ``vals`` once, accumulating the sum and count,
    then divides to produce the mean.

    Args:
        vals: An iterable of numeric values.  Must be non-empty.

    Returns:
        The arithmetic mean as a float.

    Raises:
        ZeroDivisionError: If ``vals`` is empty.
    """
    n = 0
    s = 0.0
    for i in vals:
        s += i
        n += 1
    return s / float(n)

def median(vals):
    """Compute the median of a list of numbers.

    Sorts ``vals`` and returns the middle value for odd-length lists, or
    the average of the two middle values for even-length lists.

    Args:
        vals: A sequence of numeric values.  Must be non-empty.

    Returns:
        The median value as a float.
    """
    lenvals = len(vals)
    sortvals = sorted(vals)

    if lenvals % 2 == 0:
        return (sortvals[lenvals // 2] + sortvals[lenvals // 2 - 1]) / 2.0
    else:
        return sortvals[lenvals // 2]

def mode(vals):
    """Compute the mode (most frequently occurring value) of a sequence.

    Uses :class:`collections.Counter` to count occurrences and returns
    the value with the highest count.  If multiple values share the
    maximum count, the one encountered first during dict iteration is
    returned (which is insertion order in Python 3.7+).

    Args:
        vals: An iterable of hashable values.

    Returns:
        The most frequently occurring element in ``vals``, or ``None``
        if ``vals`` is empty.
    """
    top = 0
    topkey = None
    for key, val in Counter(vals).items():
        if val > top:
            top = val
            topkey = key
    return topkey


def msqerr(vals1, vals2):
    """Compute the mean squared error between two equal-length sequences.

    Calculates the average of the squared element-wise differences::

        MSE = mean((vals1[i] - vals2[i])^2  for all i)

    Args:
        vals1: A sequence of numeric values.
        vals2: A sequence of numeric values of the same length as ``vals1``.

    Returns:
        The mean squared error as a float.

    Raises:
        AssertionError: If ``vals1`` and ``vals2`` have different lengths.
    """
    assert len(vals1) == len(vals2), "lists are not the same length"


    return mean([(vals1[i] - vals2[i]) ** 2
                 for i in range(len(vals1))])



def variance(vals):
    """Compute the sample variance of a sequence of numbers.

    Uses Bessel's correction (divides by ``n - 1``) to produce an
    unbiased estimate of the population variance::

        s^2 = sum((x - mean)^2) / (n - 1)

    Args:
        vals: A sequence of at least two numeric values.

    Returns:
        The sample variance as a float.

    Raises:
        ZeroDivisionError: If ``vals`` has fewer than 2 elements.
    """
    u = mean(vals)
    return sum((x - u)**2 for x in vals) / float(len(vals)-1)

def var(vals):
    """Alias for :func:`variance`.

    Args:
        vals: A sequence of at least two numeric values.

    Returns:
        The sample variance as a float.
    """
    return variance(vals)

def sdev(vals):
    """Compute the sample standard deviation of a sequence of numbers.

    Returns the square root of the sample variance computed by
    :func:`variance` (Bessel-corrected, ``n - 1`` denominator).

    Args:
        vals: A sequence of at least two numeric values.

    Returns:
        The sample standard deviation as a float.
    """
    return sqrt(variance(vals))

def serror(vals):
    """Compute the standard error of the mean of a sequence of numbers.

    Divides the sample standard deviation by the square root of the
    sample size::

        SE = sdev(vals) / sqrt(n)

    Args:
        vals: A sequence of at least two numeric values.

    Returns:
        The standard error of the mean as a float.
    """
    return sdev(vals) / sqrt(len(vals))

def covariance(lst1, lst2):
    """Compute the sample covariance between two equal-length sequences.

    Uses Bessel's correction (divides by ``n - 1``)::

        cov(X, Y) = sum((x - mean_x) * (y - mean_y)) / (n - 1)

    Args:
        lst1: A sequence of numeric values.
        lst2: A sequence of numeric values of the same length as ``lst1``.

    Returns:
        The sample covariance as a float.
    """
    m1 = mean(lst1)
    m2 = mean(lst2)
    tot = 0.0
    for i in range(len(lst1)):
        tot += (lst1[i] - m1) * (lst2[i] - m2)
    return tot / (len(lst1)-1)


def covmatrix(mat):
    """Compute the full pairwise sample covariance matrix for a list of sequences.

    Evaluates :func:`covariance` for every pair ``(i, j)`` of rows in
    ``mat`` (including self-covariances on the diagonal, which equal the
    sample variance of that row).

    Args:
        mat: A list of ``n`` equal-length numeric sequences (rows).

    Returns:
        A ``(n, n)`` NumPy array where element ``[i, j]`` is the sample
        covariance between ``mat[i]`` and ``mat[j]``.
    """
    size = len(mat)

    flat = [covariance(mat[i], mat[j]) for i,j in ((i,j) for i in range(size) for j in range(size))]
    return np.array(flat).reshape(size, size)

def corrmatrix(mat):
    """Compute the full pairwise Pearson correlation matrix for a list of sequences.

    Evaluates :func:`corr` for every pair ``(i, j)`` of rows in
    ``mat`` (including self-correlations of 1.0 on the diagonal).

    Args:
        mat: A list of ``n`` equal-length numeric sequences (rows).

    Returns:
        A ``(n, n)`` NumPy array where element ``[i, j]`` is the Pearson
        correlation coefficient between ``mat[i]`` and ``mat[j]``.
    """
    size = len(mat)

    flat = [corr(mat[i], mat[j]) for i,j in ((i,j) for i in range(size) for j in range(size))]
    return np.array(flat).reshape(size, size)


def corr(lst1, lst2):
    """Compute the Pearson correlation coefficient between two sequences.

    Calculates::

        r = cov(lst1, lst2) / (sdev(lst1) * sdev(lst2))

    If the denominator is zero (one or both sequences have zero variance),
    returns ``1e1000`` (effectively infinity) as a sentinel value.

    Args:
        lst1: A sequence of numeric values.
        lst2: A sequence of numeric values of the same length as ``lst1``.

    Returns:
        The Pearson correlation coefficient as a float in [-1, 1], or
        ``1e1000`` if either sequence has zero standard deviation.
    """
    num = covariance(lst1, lst2)
    denom = float(sdev(lst1) * sdev(lst2))
    if denom != 0:
        return num / denom
    else:
        return 1e1000


def qqnorm(data, plot=None):
    """Generate data for a normal quantile-quantile (Q-Q) plot.

    Sorts ``data`` and generates an equal-length sample from the standard
    normal distribution (mean 0, sigma 1), also sorted.  The two sorted
    sequences can be plotted against each other to assess normality.

    Args:
        data: A sequence of numeric values to compare against the normal
            distribution.
        plot: An optional plot object with a ``plot(x, y)`` method.  If
            provided, the Q-Q data are passed to ``plot.plot`` and the
            plot object is returned.  Defaults to ``None``.

    Returns:
        If ``plot`` is ``None``: a 2-tuple ``(data2, norm)`` where
        ``data2`` is the sorted input data and ``norm`` is a sorted
        sample from N(0, 1) of the same length.
        If ``plot`` is provided: the ``plot`` object after calling
        ``plot.plot(data2, norm)``.
    """
    data2 = sorted(data)
    norm = [random.normalvariate(0, 1) for x in range(len(data2))]
    norm.sort()

    if plot == None:
        # plotting removed (no gnuplot); return data instead
        return data2, norm
    else:
        plot.plot(data2, norm)
        return plot



def fitLine(xlist, ylist):
    """Fit a least-squares line to 2-D data and return slope and intercept.

    Uses the ordinary least-squares closed-form formula::

        slope = (sum(x*y) - n*mean_x*mean_y) / (sum(x^2) - n*mean_x^2)
        inter = mean_y - slope * mean_x

    If the denominator is zero (all x values are identical), slope is set
    to ``1e10`` as a sentinel for a vertical line.

    Args:
        xlist: A sequence of x-coordinates (numeric).
        ylist: A sequence of y-coordinates (numeric) of the same length
            as ``xlist``.

    Returns:
        A 2-tuple ``(slope, inter)`` where ``slope`` is the gradient and
        ``inter`` is the y-intercept of the fitted line.
    """
    xysum = 0
    xxsum = 0
    n = len(xlist)
    for i in range(n):
        xysum += xlist[i] * ylist[i]
        xxsum += xlist[i] * xlist[i]
    avgx = mean(xlist)
    avgy = mean(ylist)

    if (xxsum - n*avgx*avgx) == 0:
        slope = 1e10
    else:
        slope = (xysum - n*avgx*avgy) / float(xxsum - n*avgx*avgx)

    inter = (avgy*xxsum - avgx*xysum) / float(xxsum - n*avgx*avgx)

    return (slope, inter)


def fitLineError(xlist, ylist, slope, inter):
    """Compute the mean squared error of a linear fit against data.

    Evaluates the fitted line ``y_hat = slope * x + inter`` at each x
    and averages the squared residuals::

        MSE = sum((slope*x_i + inter - y_i)^2) / n

    Args:
        xlist: A sequence of x-coordinates (numeric).
        ylist: A sequence of observed y-coordinates of the same length
            as ``xlist``.
        slope: The slope of the fitted line.
        inter: The y-intercept of the fitted line.

    Returns:
        The mean squared error of the linear fit as a float.
    """
    error = 0
    n = len(xlist)

    for i in range(n):
        error += ((xlist[i]*slope + inter) - ylist[i]) ** 2
    return error / n


def pearsonsRegression(observed, expected):
    """Compute the Pearson coefficient of determination (R^2).

    Measures how well ``expected`` values explain the variance of
    ``observed``::

        R^2 = 1 - ESS / TSS

    where ``ESS = sum((observed - expected)^2)`` is the error sum of
    squares and ``TSS = sum((observed - mean(observed))^2)`` is the
    total sum of squares.

    Args:
        observed: A sequence of observed (actual) numeric values.
        expected: A sequence of predicted values of the same length as
            ``observed``.

    Returns:
        R^2 as a float.  A value of 1.0 indicates a perfect fit;
        values near 0 indicate no explanatory power; negative values
        indicate the model is worse than predicting the mean.
    """
    # error sum of squares
    ess = sum((a - b)**2 for a, b in zip(observed, expected))

    # total sum of squares
    u = mean(observed)
    tss = sum((a - u)**2 for a in observed)

    r2 = 1 - ess / tss
    return r2


def pearsonsRegressionLine(x, y, m, b):
    """Compute R^2 for data against a linear model y = m*x + b.

    Generates expected values from the line ``y = m*x + b`` and
    delegates to :func:`pearsonsRegression`.

    Args:
        x: A sequence of x-coordinates (numeric).
        y: A sequence of observed y-coordinates of the same length as
            ``x``.
        m: The slope of the reference line.
        b: The y-intercept of the reference line.

    Returns:
        R^2 as a float indicating goodness of fit of the linear model
        to the observed data.
    """
    observed = y
    expected = [m*i + b for i in x]
    return pearsonsRegression(observed, expected)



def percentile(vals, perc, rounding=-1, sort=True):
    """Return the value at a given percentile of a sequence.

    Optionally sorts ``vals`` and returns the element at index
    ``int(perc * n)`` (round down) or ``ceil(perc * n)`` (round up),
    clamped to valid list indices.

    Args:
        vals: A sequence of numeric values.
        perc: The desired percentile as a fraction in [0, 1] (e.g. 0.5
            for the median, 0.95 for the 95th percentile).
        rounding: Controls how the fractional index is resolved.
            Use ``-1`` to floor (default) or ``1`` to ceiling.
        sort: If ``True`` (default), sort ``vals`` before indexing.
            Pass ``False`` if ``vals`` is already sorted to save time.

    Returns:
        The value in ``vals`` at the requested percentile.

    Raises:
        Exception: If ``rounding`` is not ``-1`` or ``1``.
    """
    if sort:
        vals2 = sorted(vals)
    else:
        vals2 = vals
    n = len(vals2)
    if rounding == -1:
        return vals2[max(0, min(n-1, int(perc * n)))]
    elif rounding == 1:
        return vals2[max(0, min(n-1, int(ceil(perc * n))))]
    else:
        raise Exception("rounding must be 1 or -1")


def logadd(lna, lnb):
    """Add two numbers represented in log space without underflow.

    Computes ``log(exp(lna) + exp(lnb))`` in a numerically stable way::

        logadd(lna, lnb) = log(exp(lna - lnb) + 1) + lnb

    When ``lna - lnb >= 500`` the second term is negligible and ``lna``
    is returned directly to avoid overflow.

    Args:
        lna: The natural log of the first value.
        lnb: The natural log of the second value.

    Returns:
        The natural log of the sum ``exp(lna) + exp(lnb)`` as a float.
    """
    diff = lna - lnb
    if diff < 500:
        return log(exp(diff) + 1.0) + lnb
    else:
        return lna



def smooth(vals, radius):
    """Smooth a sequence by replacing each value with a local window average.

    For each position ``i``, computes the mean of the sub-list
    ``vals[i - r : i + r + 1]`` where ``r = min(i, vlen - i - 1, radius)``
    ensures the window stays within array bounds.  Values near the
    edges therefore use a smaller effective radius.

    Note:
        Not implemented as fast as possible.
        Runtime is O(len(vals) * radius).

    Args:
        vals: A sequence of numeric values.
        radius: The maximum half-width of the averaging window (the
            window spans at most ``2*radius + 1`` elements).

    Returns:
        A list of smoothed values of the same length as ``vals``.
    """
    vals2 = []
    vlen = len(vals)

    for i in range(vlen):
        radius2 = min(i, vlen - i - 1, radius)
        vals2.append(mean(vals[i-radius2:i+radius2+1]))

    return vals2




def iter_window_index(x, xdist, esp=None):
    """Iterate sliding-window index ranges over a sorted value sequence.

    Advances a window of fixed width ``xdist`` along the value axis of
    a sorted sequence ``x``, yielding the array-index bounds and value
    bounds of the window each time a point enters or exits it.

    The window boundaries are updated one step at a time: the lower
    bound advances whenever the leading point would be expelled, and the
    upper bound advances to admit the next point.

    Args:
        x: A sorted (ascending) list of numeric values.
        xdist: The width of the sliding window in the same units as
            values in ``x``.
        esp: Unused parameter retained for API compatibility.

    Yields:
        4-tuples ``(lowi, highi, low, high)`` where ``lowi`` and
        ``highi`` are the inclusive index bounds of the current window
        in ``x``, and ``low`` / ``high`` are the corresponding value
        boundaries.
    """

    vlen = len(x)
    #if esp is None:
    #    esp = min(x[i+1] - x[i] for i in range(vlen-1)
    #              if x[i+1] - x[i] > 0) / 2.0

    # simple case
    if vlen == 0:
        return

    start = x[0]
    end = x[-1]
    window = [0]

    low = start
    high = start + xdist
    lowi = 0 # inclusive
    highi = 0 # inclusive

    # move up high boundary
    while highi+1 < vlen and x[highi+1] < high:
        highi += 1

    yield (lowi, highi, low, high)

    while highi+1 < vlen:
        low_step = x[lowi] - low    # dist until expell
        high_step = x[highi+1] - high # dist until include

        # advance though duplicates
        if low_step == 0:
            lowi += 1
            continue

        if high_step == 0:
            highi += 1
            continue

        # detrmine new low high boundary
        if low_step <= high_step:
            low = x[lowi] #+ min(esp, (high_step - low_step) / 2.0)
            high = low + xdist
            lowi += 1

        if high_step <= low_step:
            highi += 1
            if highi >= vlen: break
            high = x[highi] #+ min(esp, (low_step - high_step) / 2.0)
            low = high - xdist

        assert abs((high - low) - xdist) < .001, (low, high)

        yield (lowi, highi, low, high)


def iter_window_index_step(x, size, step, minsize=0):
    """Iterate fixed-step sliding-window index ranges over a sorted value sequence.

    Advances a window of fixed width ``size`` in increments of ``step``
    along the value axis, yielding index and value bounds for each
    window position that contains at least ``minsize`` points.

    Args:
        x: A sorted (ascending) list of numeric values.
        size: The width of each window in the same units as ``x``.
        step: The distance to advance the window centre between successive
            yields.
        minsize: Minimum number of points that must be inside the window
            for it to be yielded.  Defaults to 0.

    Yields:
        4-tuples ``(lowi, highi, low, high)`` where ``lowi`` and
        ``highi`` are the inclusive index bounds of the current window
        in ``x``, and ``low`` / ``high`` are the value boundaries.
    """
    vlen = len(x)
    start = x[0]
    end = x[-1]

    low = start
    high = start + size
    i = 1

    lowi = 0
    highi = 0

    # move up high boundary
    while highi+1 < vlen and x[highi+1] < high:
        highi += 1

    while highi < vlen and high < end:
        if highi - lowi >= minsize:
            yield lowi, highi, low, high
        low = start + i * step
        high = low + size
        i += 1

        # move up low boundary
        while lowi < vlen and x[lowi] < low:
            lowi += 1

        # move up high boundary
        while highi+1 < vlen and x[highi+1] < high:
            highi += 1



def iter_window(x, xdist, func=lambda win: win, minsize=0):
    """Apply a function to each sliding window over a sorted sequence.

    Wraps :func:`iter_window_index` and yields the window midpoint
    together with ``func`` applied to the window slice.

    Note:
        The internal call uses ``xsize`` rather than ``xdist``; this is
        a latent bug in the original code and is preserved here.

    Args:
        x: A sorted (ascending) list of numeric values.
        xdist: The width of the sliding window.
        func: A callable applied to each window slice ``x[lowi:highi]``.
            Defaults to the identity function.
        minsize: Minimum number of points in the window before it is
            yielded.  Defaults to 0.

    Yields:
        2-tuples ``(midpoint, func(window))`` where ``midpoint`` is
        ``(low + high) / 2`` and ``window`` is the slice of ``x``
        within the current bounds.
    """
    for lowi, highi, low, high in iter_window_index(x, xsize):
        if highi - lowi >= minsize:
            yield (high + low)/2.0, func(x[lowi:highi])


def iter_window_step(x, width, step, func=lambda win: win, minsize=0):
    """Apply a function to each fixed-step sliding window over a sorted sequence.

    Wraps :func:`iter_window_index_step` and yields the window midpoint
    together with ``func`` applied to the window slice.  ``x`` must be
    sorted in ascending order.

    Args:
        x: A sorted (ascending) list of numeric values.
        width: The width of each window in the same units as ``x``.
        step: The distance to advance the window between successive yields.
        func: A callable applied to each window slice ``x[lowi:highi]``.
            Defaults to the identity function.
        minsize: Minimum number of points that must be in the window for
            it to be yielded.  Defaults to 0.

    Yields:
        2-tuples ``(midpoint, func(window))`` where ``midpoint`` is
        ``(low + high) / 2.0`` and ``window`` is the slice of ``x``
        within the current bounds.
    """
    for lowi, highi, low, high in iter_window_index_step(x, width, step, minsize):
        yield (high + low) / 2.0, func(x[lowi:highi])


def _sortTogether(x, y):
    """Sort two sequences together by the values of ``x``.

    Zips ``x`` and ``y`` into pairs, sorts by the first element of each
    pair, then unzips back into two separate lists.

    Args:
        x: A sequence of sortable values used as the sort key.
        y: A sequence of values of the same length as ``x``.

    Returns:
        A 2-tuple ``(x2, y2)`` where both lists have been reordered so
        that ``x2`` is sorted ascending.  Returns ``([], [])`` if ``x``
        is empty.
    """
    if not x:
        return [], []
    pairs = sorted(zip(x, y))
    x2, y2 = zip(*pairs)
    return list(x2), list(y2)


def smooth2(x, y, xradius, minsize=0, sort=False):
    """Smooth paired (x, y) data by averaging within a sliding x-radius window.

    For each point ``x[i]``, the window spans all points whose x-value
    lies within ``[x[i] - r, x[i] + r]`` where
    ``r = min(x[i] - min(x), max(x) - x[i], xradius)`` so that the
    effective radius shrinks near the data boundaries.

    Args:
        x: A sorted (ascending) list of x-coordinates.  Must be
            non-empty and of the same length as ``y``.
        y: A list of y-values corresponding to ``x``.
        xradius: The maximum half-width of the averaging window in
            the same units as ``x``.
        minsize: Minimum number of points that must be in the window
            for the averaged point to be included in the output.
            Defaults to 0.
        sort: If ``True``, sort ``x`` and ``y`` together by ``x`` before
            smoothing.  Defaults to ``False``.

    Returns:
        A 2-tuple ``(x2, y2)`` of lists containing the smoothed x and y
        values.  Returns ``([], [])`` if ``x`` is empty.
    """

    vlen = len(x)
    assert vlen == len(y)

    # simple case
    if vlen == 0:
        return [], []

    if sort:
        x, y = _sortTogether(x, y)

    x2 = []
    y2 = []

    start = min(x)
    end = max(x)
    xtot = x[0]
    ytot = y[0]

    low = 0
    high = 0

    for i in range(vlen):
        xi = x[i]

        xradius2 = min(xi - start, end - xi, xradius)

        # move window
        while x[low] < xi - xradius2:
            xtot -= x[low]
            ytot -= y[low]
            low += 1
        while x[high] < xi + xradius2:
            high += 1
            xtot += x[high]
            ytot += y[high]

        denom = float(high - low + 1)
        if denom >= minsize:
            x2.append(xtot / denom)
            y2.append(ytot / denom)

    return x2, y2


def factorial(x, k=1):
    """Compute the partial factorial product x! / k!.

    Calculates the product of all integers from ``k+1`` to ``x``
    inclusive.  When ``k=1`` (the default) this is the standard
    factorial ``x!``.  When ``k > 1`` it returns the falling factorial
    ``x! / k!``.

    Args:
        x: The upper bound of the product (inclusive).  Converted to
            ``int`` internally.
        k: The lower bound; the product starts at ``k+1``.  Defaults
            to 1.

    Returns:
        An integer equal to ``(k+1) * (k+2) * ... * x``, or 1 if the
        range is empty (i.e. ``x <= k``).
    """
    n = 1
    for i in range(int(k)+1, int(x)+1):
        n *= i
    return n


def logfactorial(x, k=1):
    """Compute log(x! / k!) in log space.

    Returns the natural log of the partial factorial product
    ``(k+1) * (k+2) * ... * x`` by summing ``log(i)`` terms.  This
    avoids integer overflow for large ``x``.

    Args:
        x: The upper bound of the product (inclusive).  Converted to
            ``int`` internally.
        k: The lower bound; the product starts at ``k+1``.  Defaults
            to 1.

    Returns:
        A float equal to ``log((k+1)) + log(k+2) + ... + log(x)``,
        or 0.0 if the range is empty.
    """
    n = 0
    for i in range(int(k)+1, int(x)+1):
        n += log(i)
    return n


def choose(n, k):
    """Compute the binomial coefficient C(n, k) = n! / (k! * (n-k)!).

    Uses a multiplicative formula for efficiency, exploiting the
    symmetry ``C(n, k) == C(n, n-k)`` to minimise the number of
    multiplications.  Returns the result rounded to the nearest integer.

    Args:
        n: The total number of items.
        k: The number of items to choose.

    Returns:
        An integer equal to C(n, k).  Returns 1.0 when both ``n`` and
        ``k`` are 0, and 0 when any argument is negative or ``k > n``.
    """
    if n == 0 and k == 0:
        return 1.0

    if n < 0 or k < 0 or k > n:
        return 0

    # optimization for speed
    if k > n/2:
        k = n - k

    t = 1.0
    for i in range(1, k+1):
        t = t * (n - i + 1) / i
    return int(t + 0.5)
    #return factorial(n, n - k) / factorial(k)


def _oneNorm(weights):
    """Normalise a list of weights so they sum to 1.

    Divides each weight by the total sum of all weights.

    Args:
        weights: A list of non-negative numeric values whose sum is
            positive.

    Returns:
        A new list of floats of the same length as ``weights`` that
        sum to 1.0.
    """
    s = sum(weights)
    return [w / s for w in weights]


def sample(weights):
    """Randomly choose an index proportional to the given weights.

    Normalises ``weights`` to a proper probability distribution and then
    samples using a CDF built from the normalised weights and a binary
    search via :func:`algorithms.binsearch`.

    Item ``i`` is chosen with probability ``weights[i] / sum(weights)``.

    Args:
        weights: A list of non-negative numeric values.  The length
            determines the range of possible return values (0 to
            ``len(weights) - 1``).

    Returns:
        An integer index into ``weights``, selected with probability
        proportional to each weight.

    Raises:
        AssertionError: If ``algorithms.binsearch`` returns ``None`` for
            the lower bound, indicating an unexpected state.
    """
    probs = _oneNorm(weights)

    cdf = [0]
    for i in range(1, len(probs)):
        cdf.append(cdf[-1] + probs[i-1])

    pick = random.random()

    low,top = algorithms.binsearch(cdf, pick)

    assert low != None

    return low


def chyper(m, n, M, N, report=0):
    """Compute a hypergeometric cumulative probability via an external ``chyper`` binary.

    Models drawing ``n`` balls from an urn containing ``N`` balls of which
    ``M`` are white (successes).  ``m`` is the number of white balls drawn.
    Calls the external command-line tool ``chyper`` and parses its output.

    Args:
        m: Number of white balls drawn (observed successes).  Must be
            an ``int`` with ``m <= n`` and ``m <= M``.
        n: Total balls drawn.  Must be an ``int`` with ``n <= N``.
        M: Total white balls in urn.  Must be an ``int``.
        N: Total balls in urn.  Must be an ``int``.
        report: Controls which tail(s) are returned.
            ``0`` — p-value for over-representation (default).
            ``1`` — p-value for under-representation.
            ``2`` — 2-tuple ``(over_p, under_p)``.

    Returns:
        A float p-value, or a list of two floats when ``report=2``.

    Raises:
        AssertionError: If arguments do not satisfy type or range constraints.
        Exception: If the ``chyper`` command produces no output.
        Exception: If ``report`` is not 0, 1, or 2.
    """

    assert( (type(m) == type(n) == type(M) == type(N) == int)
            and m <= n and m <= M and n <= N)

    command = "chyper %d %d %d %d 2>/dev/null" % (m, n, M, N)
    stream = os.popen(command)
    val = stream.read()
    if val == '':
        raise Exception("error in chyper")
    else:
        val = val.strip()
        vals = list(map(float, val.split(' ')[4:6]))

    if report == 0:
        #p-val for over-repr.
        return vals[0]
    elif report == 1:
        #p-val for under-repr.
        return vals[1]
    elif report == 2:
        #tuple (over, under)
        return vals
    else:
        raise Exception("unknown option")


def rhyper(m, n, M, N, report=0):
    """Compute a hypergeometric cumulative probability via R (rpy2).

    Models drawing ``n`` balls from an urn containing ``N`` balls of which
    ``M`` are white (successes).  ``m`` is the number of white balls drawn.
    Uses R's ``phyper`` function via the rpy2 interface.

    Args:
        m: Number of white balls drawn (observed successes).  Must be
            an ``int`` with ``m <= n`` and ``m <= M``.
        n: Total balls drawn.  Must be an ``int`` with ``n <= N``.
        M: Total white balls in urn.  Must be an ``int``.
        N: Total balls in urn.  Must be an ``int``.
        report: Controls which tail(s) are returned.
            ``0`` — p-value for over-representation, i.e.
            ``P(X >= m)`` (default).
            ``1`` — p-value for under-representation, i.e. ``P(X <= m)``.
            ``2`` — 2-tuple ``(over_p, under_p)``.

    Returns:
        A float p-value, or a 2-tuple of floats when ``report=2``.

    Raises:
        AssertionError: If arguments do not satisfy type or range constraints.
        Exception: If ``report`` is not 0, 1, or 2.
    """

    import rpy2.robjects as r_module
    r = r_module.r

    assert( (type(m) == type(n) == type(M) == type(N) == int)
            and m <= n and m <= M and n <= N)

    if report == 0:
        #p-val for over-repr.
        return r['phyper'](m-1, M, N-M, n, **{'lower.tail': False})[0]
    elif report == 1:
        #p-val for under-repr.
        return r['phyper'](m, M, N-M, n)[0]
    elif report == 2:
        #tuple (over, under)
        return r['phyper'](m-1, M, N-M, n, **{'lower.tail': False})[0], r['phyper'](m, M, N-M, n)[0]
    else:
        raise Exception("unknown option")

def cdf(vals):
    """Compute the empirical cumulative distribution function (ECDF) of a list.

    Sorts ``vals`` and assigns each unique value a cumulative probability
    equal to its 0-based rank divided by the total number of values.

    Args:
        vals: A sequence of numeric values.

    Returns:
        A 2-tuple ``(x, y)`` where ``x`` is the sorted list of values and
        ``y`` is the corresponding list of cumulative probabilities in
        [0, 1).
    """
    vals = sorted(vals)
    tot = float(len(vals))
    x = []
    y = []

    for i, x2 in enumerate(vals):
        x.append(x2)
        y.append(i / tot)

    return x, y


def enrichItems(in_items, out_items, M=None, N=None, useq=True, extra=False):
    """Calculate item enrichment between an in-set and an out-set.

    Counts how often each item appears in ``in_items`` vs ``out_items`` and
    tests for enrichment using the hypergeometric distribution via
    :func:`rhyper`.  Optionally adjusts p-values to q-values (FDR) and
    adds fold-enrichment columns.

    Args:
        in_items: An iterable of items in the foreground (in-set).
        out_items: An iterable of items in the background (out-set).
        M: The foreground population size.  Defaults to
            ``len(in_items)``.
        N: The total population size.  Defaults to
            ``len(in_items) + len(out_items)``.
        useq: If ``True`` (default), add ``qval`` and ``qval_under``
            columns computed via FDR correction using :func:`qvalues`.
        extra: If ``True``, add columns ``in_size``, ``out_size``,
            ``item_ratio``, ``size_ratio``, and ``fold`` for fold-
            enrichment analysis.  Defaults to ``False``.

    Returns:
        A :class:`pandas.DataFrame` sorted by ``pval`` (ascending) with
        columns ``item``, ``in_count``, ``out_count``, ``pval``,
        ``pval_under``, and optionally ``qval``, ``qval_under``, and
        fold-enrichment columns.
    """

    # count items using defaultdict instead of rasmus util.Dict
    counts = defaultdict(lambda: [0, 0])
    for item in in_items:
        counts[item][0] += 1
    for item in out_items:
        counts[item][1] += 1

    if N is None:
        N = len(in_items) + len(out_items)
    if M is None:
        M = len(in_items)

    rows = []
    for item, (a, b) in counts.items():
        rows.append(dict(
            item=item,
            in_count=a,
            out_count=b,
            pval=rhyper(a, a+b, M, N),
            pval_under=rhyper(a, a+b, M, N, 1)
        ))

    tab = pd.DataFrame(rows, columns=["item", "in_count", "out_count", "pval", "pval_under"])

    # add qvalues
    if useq:
        qval = qvalues(list(tab["pval"]))
        qval_under = qvalues(list(tab["pval_under"]))

        tab["qval"] = qval
        tab["qval_under"] = qval_under

    if extra:
        tab["in_size"] = M
        tab["out_size"] = N - M
        tab["item_ratio"] = tab.apply(
            lambda row: row["in_count"] / float(row["in_count"] + row["out_count"]), axis=1)
        tab["size_ratio"] = M / float(N)
        tab["fold"] = tab["item_ratio"] / tab["size_ratio"]

    tab = tab.sort_values("pval").reset_index(drop=True)
    return tab


def qvalues(pvals):
    """Compute Benjamini-Hochberg FDR-adjusted p-values (q-values) via R.

    Calls R's ``p.adjust`` function with ``method='fdr'`` through rpy2.

    Args:
        pvals: A list of raw p-values (floats in [0, 1]).

    Returns:
        A list of FDR-adjusted p-values (q-values) of the same length
        as ``pvals``.
    """
    import rpy2.robjects as robjects
    ret = robjects.r['p.adjust'](robjects.FloatVector(pvals), 'fdr')
    return list(ret)

def qvalues2(pvals):
    """Compute q-values using the Storey-Tibshirani method via R's qvalue package.

    Loads the ``qvalue`` R package through rpy2 and calls ``qvalue()`` on
    the provided p-values.

    Args:
        pvals: A list of raw p-values (floats in [0, 1]).

    Returns:
        A list of q-values of the same length as ``pvals`` as computed
        by the Storey-Tibshirani estimator.
    """
    import rpy2.robjects as robjects
    robjects.r['library']('qvalue')
    ret = robjects.r['qvalue'](robjects.FloatVector(pvals))
    return list(ret.rx2('qvalues'))


#=============================================================================
# Distributions
#

def uniformPdf(x, params):
    """Evaluate the Uniform(a, b) probability density function at ``x``.

    Returns ``1 / (b - a)`` when ``a <= x <= b``, and 0 otherwise.

    Args:
        x: The point at which to evaluate the PDF.
        params: A 2-tuple ``(a, b)`` defining the lower and upper bounds
            of the uniform distribution.

    Returns:
        The PDF value at ``x`` as a float.
    """
    a, b = params
    if x < a or x > b:
        return 0.0
    else:
        return 1.0 / (b - a)


def binomialPdf(k, params):
    """Evaluate the Binomial(n, p) probability mass function at ``k``.

    Computes::

        P(X = k) = C(n, k) * p^k * (1 - p)^(n - k)

    Args:
        k: The number of successes (non-negative integer).
        params: A 2-tuple ``(p, n)`` where ``p`` is the success probability
            per trial and ``n`` is the total number of trials.

    Returns:
        The probability of exactly ``k`` successes as a float.
    """
    p, n = params
    return choose(n, k) * (p ** k) * ((1.0-p) ** (n - k))

def gaussianPdf(x, params):
    """Evaluate the standard Normal N(0, 1) probability density function at ``x``.

    Computes::

        f(x) = (1 / sqrt(2*pi)) * exp(-x^2 / 2)

    Note:
        The ``params`` argument is accepted but ignored; this function
        always evaluates the standard normal (mean 0, variance 1).

    Args:
        x: The point at which to evaluate the PDF.
        params: Unused.  Accepted for API consistency with other PDF
            functions.

    Returns:
        The standard normal PDF value at ``x`` as a float.
    """
    return 1/sqrt(2*pi) * exp(- x**2 / 2.0)

def normalPdf(x, params):
    """Evaluate the Normal(mu, sigma) probability density function at ``x``.

    Computes::

        f(x) = (1 / (sigma * sqrt(2*pi))) * exp(-(x - mu)^2 / (2*sigma^2))

    Args:
        x: The point at which to evaluate the PDF.
        params: A 2-tuple ``(mu, sigma)`` — the mean and standard
            deviation of the normal distribution.

    Returns:
        The normal PDF value at ``x`` as a float.
    """
    mu, sigma = params
    return 1.0/(sigma * sqrt(2.0*pi)) * exp(- (x - mu)**2 / (2.0 * sigma**2))

def normalCdf(x, params):
    """Evaluate the Normal(mu, sigma) cumulative distribution function at ``x``.

    Computes::

        F(x) = (1 + erf((x - mu) / (sigma * sqrt(2)))) / 2

    Args:
        x: The point at which to evaluate the CDF.
        params: A 2-tuple ``(mu, sigma)`` — the mean and standard
            deviation of the normal distribution.

    Returns:
        The cumulative probability P(X <= x) as a float in [0, 1].
    """
    mu, sigma = params
    return (1 + erf((x - mu)/(sigma * sqrt(2)))) / 2.0

def logNormalPdf(x, params):
    """Evaluate the log-normal probability density function at ``x``.

    The log-normal distribution describes a variable whose natural
    logarithm is normally distributed.  The PDF is::

        f(x) = (1 / (x * sigma * sqrt(2*pi))) * exp(-(log(x) - mu)^2 / (2*sigma^2))

    Args:
        x: The point at which to evaluate the PDF.  Must be positive.
        params: A 2-tuple ``(mu, sigma)`` — the mean and standard
            deviation of the variable's natural logarithm.

    Returns:
        The log-normal PDF value at ``x`` as a float.  Returns nonsensical
        values for ``x <= 0``.
    """
    mu, sigma = params
    return 1/(x * sigma * sqrt(2*pi)) * \
           exp(- (log(x) - mu)**2 / (2.0 * sigma**2))

def logNormalCdf(x, params):
    """Evaluate the log-normal cumulative distribution function at ``x``.

    Computes::

        F(x) = (1 + erf((log(x) - mu) / (sigma * sqrt(2)))) / 2

    Args:
        x: The point at which to evaluate the CDF.  Must be positive.
        params: A 2-tuple ``(mu, sigma)`` — the mean and standard
            deviation of the variable's natural logarithm.

    Returns:
        The cumulative probability P(X <= x) as a float in [0, 1].
    """
    mu, sigma = params
    return (1 + erf((log(x) - mu)/(sigma * sqrt(2)))) / 2.0


def poissonPdf(x, params):
    """Evaluate the Poisson probability mass function at ``x``.

    Computes the probability in log space to avoid overflow::

        P(X = x) = exp(-lambda) * lambda^x / x!
                 = exp(-lambda + sum(log(lambda/i) for i in 1..x))

    Args:
        x: The number of events (non-negative integer).
        params: A 1-tuple or list whose first element is ``lambda``
            (the expected number of events, must be positive).

    Returns:
        The Poisson PMF value P(X = x) as a float.  Returns 0.0 if
        ``x < 0`` or ``lambda <= 0``.
    """
    lambd = params[0]

    if x < 0 or lambd <= 0:
        return 0.0

    a = 0
    for i in range(1, int(x)+1):
        a += log(lambd / float(i))
    return exp(-lambd + a)


def poissonCdf(x, params):
    """Evaluate the Poisson cumulative distribution function at ``x``.

    Computes P(X <= x) using the regularised incomplete gamma function::

        F(x; lambda) = (Gamma(floor(x+1)) - gammainc(floor(x+1), lambda))
                       / floor(x)!

    Note:
        Not implemented accurately for large ``x`` or ``lambda``.

    Args:
        x: The upper bound (non-negative number; floor is taken
            internally).
        params: A 1-tuple or list whose first element is ``lambda``
            (the expected number of events).

    Returns:
        The cumulative probability P(X <= x) as a float, or 0 if
        ``x < 0``.
    """
    # NOTE: not implemented accurately for large x or lambd
    lambd = params[0]

    if x < 0:
        return 0
    else:
        return (gamma(floor(x+1)) - gammainc(floor(x + 1), lambd)) / \
               factorial(floor(x))


def poissonvariate(lambd):
    """Draw a random sample from a Poisson distribution.

    Uses Knuth's algorithm: generate uniform random variables and
    multiply them together until their product falls below
    ``exp(-lambda)``.  The count of multiplications minus one is the
    Poisson variate.

    Args:
        lambd: The expected number of events per interval (lambda > 0).

    Returns:
        A non-negative integer drawn from Poisson(lambda).
    """
    l = exp(-lambd)
    k = 0
    p = 1.0

    while 1:
        k += 1
        p *= random.random()
        if p < l:
            return k - 1

def exponentialPdf(x, params):
    """Evaluate the Exponential(lambda) probability density function at ``x``.

    Computes::

        f(x; lambda) = lambda * exp(-lambda * x)   for x >= 0, lambda >= 0

    Args:
        x: The point at which to evaluate the PDF.
        params: A 1-tuple or list whose first element is ``lambda``
            (the rate parameter).

    Returns:
        The exponential PDF value at ``x`` as a float.  Returns 0.0 if
        ``x < 0`` or ``lambda < 0``.
    """
    lambd = params[0]

    if x < 0 or lambd < 0:
        return 0.0
    else:
        return lambd * exp(-lambd * x)


def exponentialCdf(x, params):
    """Evaluate the Exponential(lambda) cumulative distribution function at ``x``.

    Computes::

        F(x; lambda) = 1 - exp(-lambda * x)   for x >= 0, lambda >= 0

    Args:
        x: The point at which to evaluate the CDF.
        params: A 1-tuple or list whose first element is ``lambda``
            (the rate parameter).

    Returns:
        The cumulative probability P(X <= x) as a float.  Returns 0.0 if
        ``x < 0`` or ``lambda < 0``.
    """
    lambd = params[0]

    if x < 0 or lambd < 0:
        return 0.0
    else:
        return 1.0 - exp(-lambd * x)


def exponentialvariate(lambd):
    """Draw a random sample from an Exponential(lambda) distribution.

    Uses the inverse CDF (quantile) method: if U ~ Uniform(0,1) then
    ``-log(U) / lambda`` is Exponentially distributed with rate ``lambda``.

    Args:
        lambd: The rate parameter (lambda > 0).

    Returns:
        A non-negative float drawn from Exponential(lambda).
    """
    return -log(random.random()) / lambd

def gammaPdf(x, params):
    """Evaluate the Gamma(alpha, beta) probability density function at ``x``.

    Uses the rate (inverse-scale) parameterisation::

        f(x; alpha, beta) = beta^alpha * x^(alpha-1) * exp(-beta*x)
                            / Gamma(alpha)

    Args:
        x: The point at which to evaluate the PDF.  Must be positive.
        params: A 2-tuple ``(alpha, beta)`` — the shape and rate
            parameters.  Both must be positive.

    Returns:
        The gamma PDF value at ``x`` as a float.  Returns 0.0 if any of
        ``x``, ``alpha``, or ``beta`` is non-positive.
    """
    alpha, beta = params
    if x <= 0 or alpha <= 0 or beta <= 0:
        return 0.0
    else:
        return (exp(-x * beta) * (x ** (alpha - 1)) * (beta ** alpha)) / \
           gamma(alpha)

def gammaPdf2(x, params):
    """Evaluate the Gamma(alpha, beta) PDF at ``x`` using log-space arithmetic.

    Numerically more stable than :func:`gammaPdf` for large parameter
    values.  Computes the same distribution in log space::

        log f = -x*beta + (alpha-1)*log(x) + alpha*log(beta) - gammaln(alpha)

    Args:
        x: The point at which to evaluate the PDF.  Must be positive.
        params: A 2-tuple ``(alpha, beta)`` — the shape and rate
            parameters (rate parameterisation).  Both must be positive.

    Returns:
        The gamma PDF value at ``x`` as a float.  Returns 0.0 if any of
        ``x``, ``alpha``, or ``beta`` is non-positive.
    """
    alpha, beta = params
    if x <= 0 or alpha <= 0 or beta <= 0:
        return 0.0
    else:
        return exp(- x * beta + (alpha - 1)*log(x) + alpha * log(beta) -
                   gammaln(alpha))


def gammaCdf(x, params):
    """Evaluate the Gamma(alpha, beta) cumulative distribution function at ``x``.

    Computes P(X <= x) using the lower incomplete gamma function::

        F(x; alpha, beta) = gammainc(alpha, x*beta) / Gamma(alpha)

    Args:
        x: The point at which to evaluate the CDF.
        params: A 2-tuple ``(alpha, beta)`` — the shape and rate
            parameters (rate parameterisation).  Both must be positive.

    Returns:
        The cumulative probability P(X <= x) as a float.  Returns 0 if
        ``x <= 0``.
    """
    alpha, beta = params
    if x <= 0:
        return 0
    else:
        return gammainc(alpha, x * beta) / gamma(alpha)


def betaPdf2(x, params):
    """Evaluate the Beta(alpha, beta) PDF at ``x`` using direct gamma computation.

    Simpler but less numerically stable than :func:`betaPdf`; will
    overflow for ``alpha`` or ``beta`` values near 100 because it
    evaluates ``Gamma(alpha + beta)`` directly.

    Formula::

        f(x; alpha, beta) = Gamma(alpha+beta) / (Gamma(alpha)*Gamma(beta))
                            * x^(alpha-1) * (1-x)^(beta-1)

    Args:
        x: The point at which to evaluate the PDF.  Must satisfy
            ``0 < x < 1``.
        params: A 2-tuple ``(alpha, beta)`` — the shape parameters, both
            must be positive.

    Returns:
        The beta PDF value at ``x`` as a float.  Returns 0.0 if ``x``
        is outside (0, 1) or if either shape parameter is non-positive.
    """
    alpha, beta = params
    if 0 < x < 1 and alpha > 0 and beta > 0:
        return gamma(alpha + beta) / (gamma(alpha)*gamma(beta)) * \
               x ** (alpha-1) * (1-x)**(beta-1)
    else:
        return 0.0

def betaPdf(x, params):
    """Evaluate the Beta(alpha, beta) PDF at ``x`` using log-gamma arithmetic.

    Numerically stable implementation that avoids overflow by computing
    the PDF in log space::

        log f = gammaln(alpha+beta) - gammaln(alpha) - gammaln(beta)
                + (alpha-1)*log(x) + (beta-1)*log(1-x)

    Args:
        x: The point at which to evaluate the PDF.  Must satisfy
            ``0 < x < 1``.
        params: A 2-tuple ``(alpha, beta)`` — the shape parameters, both
            must be positive.

    Returns:
        The beta PDF value at ``x`` as a float.  Returns 0.0 if ``x``
        is outside (0, 1) or if either shape parameter is non-positive.
    """
    alpha, beta = params

    if 0 < x < 1 and alpha > 0 and beta > 0:
        return e**(gammaln(alpha + beta) - (gammaln(alpha) + gammaln(beta)) + \
                   (alpha-1) * log(x) +  (beta-1) * log(1-x))
    else:
        return 0.0



def betaPdf3(x, params):
    alpha, beta = map(int, params)
    if 0 < x < 1 and alpha > 0 and beta > 0:
        n = min(alpha-1, beta-1)
        m = max(alpha-1, beta-1)

        prod1 = 1
        for i in range(1,n+1):
            prod1 *= ((n+i)*x*(1-x))/i

        prod2 = 1
        if alpha > beta:
            for i in range(n+1, m+1):
                prod2 *= ((n+i)*x)/i
        else:
            for i in range(n+1, m+1):
                prod2 *= ((n+i)*(1-x))/i

        return prod1 * prod2 * (alpha + beta - 1)
    else:
        return 0.0


def gamma(x):
    """
    Lanczos approximation to the gamma function.

    found on http://www.rskey.org/gamma.htm
    """

    ret = 1.000000000190015 + \
          76.18009172947146 / (x + 1) + \
          -86.50532032941677 / (x + 2) + \
          24.01409824083091 / (x + 3) + \
          -1.231739572450155 / (x + 4) + \
          1.208650973866179e-3 / (x + 5) + \
          -5.395239384953e-6 / (x + 6)

    return ret * sqrt(2*pi)/x * (x + 5.5)**(x+.5) * exp(-x-5.5)



def gammaln(xx):
    """
    From numerical alogrithms in C

    float gammln(float xx)
    Returns the value ln[(xx)] for xx > 0.
    {
        Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
        accuracy is good enough.
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
             24.01409824083091,-1.231739572450155,
             0.1208650973866179e-2,-0.5395239384953e-5};
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
    }
    """

    cof = [76.18009172947146,-86.50532032941677,
         24.01409824083091,-1.231739572450155,
         0.1208650973866179e-2,-0.5395239384953e-5]

    y = x = xx
    tmp = x + 5.5
    tmp -= (x + 0.5) * log(tmp)
    ser = 1.000000000190015

    for j in range(6):
        y += 1
        ser += cof[j] / y

    return - tmp + log(2.5066282746310005 * ser / x)




GAMMA_INCOMP_ACCURACY = 1000
def gammainc(a, x):
    """Lower incomplete gamma function"""
    # found on http://www.rskey.org/gamma.htm

    ret = 0
    term = 1.0/x
    for n in range(GAMMA_INCOMP_ACCURACY):
        term *= x/(a+n)
        ret += term
        if term < .0001:
            break
    return x**a * exp(-x) * ret


def erf(x):
    # http://www.theorie.physik.uni-muenchen.de/~serge/erf-approx.pdf

    a = 8/(3*pi) * (pi - 3)/(4 - pi)
    axx = a * x * x

    if x >= 0:
        return sqrt(1 - exp(-x*x * (4.0/pi + axx)/(1 + axx)))
    else:
        return - sqrt(1 - exp(-x*x * (4.0/pi + axx)/(1 + axx)))



def chiSquare(rows, expected=None, nparams=0):
    # ex: rows = [[1,2,3],[1,4,5]]
    assert(len(set(map(len, rows))) <= 1)

    if 0 in map(sum,rows): return 0,1.0
    cols = zip(* rows)
    if 0 in map(sum,cols): return 0,1.0

    if not expected:
        expected = make_expected(rows)

    chisq = 0
    for obss,exps in zip(rows,expected):
        for obs, exp in zip(obss, exps):
            chisq += ((obs-exp)**2)/exp

    df = max(len(rows)-1, 1)*max(len(rows[0])-1, 1) - nparams

    p = chi_square_lookup(chisq,df)

    return chisq,p


def make_expected(rows):
    rowtotals = map(sum, rows)
    coltotals = map(sum, zip(* rows))
    grandtotal = float(sum(rowtotals))

    expected = []
    for row,rowtotal in zip(rows,rowtotals):
        expected_row = []
        for obs, coltotal in zip(row, coltotals):
            exp = rowtotal * coltotal / grandtotal
            expected_row.append(exp)
        expected.append(expected_row)
    return expected


def chiSquareFit(xbins, ybins, func, nsamples, nparams, minsamples=5):
    sizes = [xbins[i+1] - xbins[i] for i in range(len(xbins)-1)]
    sizes.append(sizes[-1])

    # only focus on bins that are large enough
    counts = [ybins[i] * sizes[i] * nsamples for i in range(len(xbins)-1)]

    expected = []
    for i in range(len(xbins)-1):
        expected.append((func(xbins[i]) + func(xbins[i+1]))/2.0 *
                         sizes[i] * nsamples)

    # ensure we have enough expected samples in each bin
    ind = [i for i, v in enumerate(expected) if v >= minsamples]
    counts = [counts[i] for i in ind]
    expected = [expected[i] for i in ind]

    if len(counts) == 0:
        return [0, 1], counts, expected
    else:
        return chiSquare([counts], [expected], nparams), counts, expected


chi_square_table = {
    1: [1.64, 2.71, 3.84, 5.02, 6.64, 10.83],
    2: [3.22, 4.61, 5.99, 7.38, 9.21, 13.82],
    3: [4.64, 6.25, 7.82, 9.35, 11.34, 16.27],
    4: [5.99, 7.78, 9.49, 11.14, 13.28, 18.47],
    5: [7.29, 9.24, 11.07, 12.83, 15.09, 20.52],
    6: [8.56, 10.64, 12.59, 14.45, 16.81, 22.46],
    7: [9.80, 12.02, 14.07, 16.01, 18.48, 24.32],
    8: [11.03, 13.36, 15.51, 17.53, 20.09, 26.12],
    9: [12.24, 14.68, 16.92, 19.02, 21.67, 27.88],
    10: [13.44, 15.99, 18.31, 20.48, 23.21, 29.59],
    11: [14.63, 17.28, 19.68, 21.92, 24.72, 31.26],
    12: [15.81, 18.55, 21.03, 23.34, 26.22, 32.91],
    13: [16.98, 19.81, 22.36, 24.74, 27.69, 34.53],
    14: [18.15, 21.06, 23.68, 26.12, 29.14, 36.12],
    15: [19.31, 22.31, 25.00, 27.49, 30.58, 37.70],
    16: [20.47, 23.54, 26.30, 28.85, 32.00, 39.25],
    17: [21.61, 24.77, 27.59, 30.19, 33.41, 40.79],
    18: [22.76, 25.99, 28.87, 31.53, 34.81, 42.31],
    19: [23.90, 27.20, 30.14, 32.85, 36.19, 43.82],
    20: [25.04, 28.41, 31.41, 34.17, 37.57, 45.31],
    21: [26.17, 29.62, 32.67, 35.48, 38.93, 46.80],
    22: [27.30, 30.81, 33.92, 36.78, 40.29, 48.27],
    23: [28.43, 32.01, 35.17, 38.08, 41.64, 49.73],
    24: [29.55, 33.20, 36.42, 39.36, 42.98, 51.18],
    25: [30.68, 34.38, 37.65, 40.65, 44.31, 52.62],
    26: [31.79, 35.56, 38.89, 41.92, 45.64, 54.05],
    27: [32.91, 36.74, 40.11, 43.19, 46.96, 55.48],
    28: [34.03, 37.92, 41.34, 44.46, 48.28, 56.89],
    29: [35.14, 39.09, 42.56, 45.72, 49.59, 58.30],
    30: [36.25, 40.26, 43.77, 46.98, 50.89, 59.70]
}


def chi_square_lookup(value, df):

    ps = [0.20, 0.10, 0.05, 0.025, 0.01, 0.001]

    if df <= 0:
        return 1.0

    row = chi_square_table[min(df, 30)]

    for i in range(0,len(row)):
        if row[i] >= value:
            i = i-1
            break

    if i == -1: return 1
    else: return ps[i]


def ttest(lst1, lst2):
    sdevdist = sqrt(var(lst1)/len(lst1) + var(lst2)/len(lst2))
    t = abs(mean(lst1) - mean(lst2)) / sdevdist
    df = len(lst2) + len(lst2) - 2

"""
t-table

 	0.1  	0.05  	0.01  	0.001
1 	6.31 	12.71 	63.66 	636.62
2 	2.92 	4.30 	9.93 	31.60
3 	2.35 	3.18 	5.84 	12.92
4 	2.13 	2.78 	4.60 	8.61
5 	2.02 	2.57 	4.03 	6.87
6 	1.94 	2.45 	3.71 	5.96
7 	1.89 	2.37 	3.50 	5.41
8 	1.86 	2.31 	3.36 	5.04
9 	1.83 	2.26 	3.25 	4.78
10 	1.81 	2.23 	3.17 	4.59
11 	1.80 	2.20 	3.11 	4.44
12 	1.78 	2.18 	3.06 	4.32
13 	1.77 	2.16 	3.01 	4.22
14 	1.76 	2.14 	2.98 	4.14
15 	1.75 	2.13 	2.95 	4.07
16 	1.75 	2.12 	2.92 	4.02
17 	1.74 	2.11 	2.90 	3.97
18 	1.73 	2.10 	2.88 	3.92
19 	1.73 	2.09 	2.86 	3.88
20 	1.72 	2.09 	2.85 	3.85
21 	1.72 	2.08 	2.83 	3.82
22 	1.72 	2.07 	2.82 	3.79
23 	1.71 	2.07 	2.82 	3.77
24 	1.71 	2.06 	2.80 	3.75
25 	1.71 	2.06 	2.79 	3.73
26 	1.71 	2.06 	2.78 	3.71
27 	1.70 	2.05 	2.77 	3.69
28 	1.70 	2.05 	2.76 	3.67
29 	1.70 	2.05 	2.76 	3.66
30 	1.70 	2.04 	2.75 	3.65
40 	1.68 	2.02 	2.70 	3.55
60 	1.67 	2.00 	2.66 	3.46
120 1.66 	1.98 	2.62 	3.37
"""

"""
r	90%	95%	97.5%	99.5%
1	3.07766	6.31371	12.7062	63.656
2	1.88562	2.91999	4.30265	9.92482
3	1.63774	2.35336	3.18243	5.84089
4	1.53321	2.13185	2.77644	4.60393
5	1.47588	2.01505	2.57058	4.03212
10	1.37218	1.81246	2.22814	3.16922
30	1.31042	1.69726	2.04227	2.74999
100	1.29007	1.66023	1.98397	2.62589
infty	1.28156	1.64487	1.95999	2.57584
"""


def spearman(vec1, vec2):
    """Spearman's rank test"""

    assert len(vec1) == len(vec2), "vec1 and vec2 are not the same length"

    n = len(vec1)
    rank1 = sorted(range(len(vec1)), key=lambda i: vec1[i])
    rank2 = sorted(range(len(vec2)), key=lambda i: vec2[i])

    R = sum((vec1[i] - vec2[i])**2 for i in range(n))

    Z = (6*R - n*(n*n - 1)) / (n*(n + 1) * sqrt(n - 1))

    return Z



# input:
#   xdata, ydata  - data to fit
#   func          - a function of the form f(x, params)
#
def fitCurve(xdata, ydata, func, paramsInit):
    import scipy.optimize

    y = np.array(ydata)
    p0 = np.array(paramsInit)

    def error(params):
        y2 = np.array([func(x, params) for x in xdata])
        return y - y2

    params, msg = scipy.optimize.leastsq(error, p0)

    resid = error(params)

    return list(params), sum(resid*resid)


def fitDistrib(func, paramsInit, data, start, end, step, perc=1.0):
    # NOTE: fitDistrib is disabled because it depends on rasmus util.distrib
    # and util.histbins which are not available.
    # xdata, ydata = util.distrib(data, low=start, width=step)
    # ydata = [i / perc for i in ydata]
    # xdata = util.histbins(xdata)
    # params, resid = fitCurve(xdata, ydata, func, paramsInit)
    # return params, resid
    raise NotImplementedError("fitDistrib requires rasmus util.distrib which is not available")


def plotfuncFit(func, paramsInit, xdata, ydata, start, end, step, plot=None,
                **options):
    # NOTE: plotting via gnuplot removed; returns params and resid only
    params, resid = fitCurve(xdata, ydata, func, paramsInit)
    # plot.plot(util.histbins(xdata), ydata, **options)
    # plot.plotfunc(lambda x: func(x, params), start, end, step)
    return None, params, resid


def plotdistribFit(func, paramsInit, data, start, end, step, plot=None,
                   **options):
    # NOTE: disabled because it requires rasmus util.distrib
    raise NotImplementedError("plotdistribFit requires rasmus util.distrib which is not available")



def solveCubic(a, b, c, real=True):
    """solves x^3 + ax^2 + bx + c = 0 for x"""

    p = b - a*a / 3.0
    q = c + (2*a*a*a - 9*a*b) / 27.0

    # special case: avoids division by zero later on
    if p == q == 0:
        return [- a / 3.0]

    #
    # u = (q/2 +- sqrt(q^2/4 + p^3/27))^(1/3)
    #

    # complex math is used to find complex roots
    sqrteqn = cmath.sqrt(q*q/4.0 + p*p*p/27.0)

    # find fist cube root
    u1 = (q/2.0 + sqrteqn)**(1/3.0)

    # special case: avoids division by zero later on
    if u1 == 0:
        u1 = (q/2.0 - sqrteqn)**(1/3.0)

    # find other two cube roots
    u2 = u1 * complex(-.5, -sqrt(3)/2)
    u3 = u1 * complex(-.5, sqrt(3)/2)

    # finds roots of cubic polynomial
    root1 = p / (3*u1) - u1 - a / 3.0
    root2 = p / (3*u2) - u2 - a / 3.0
    root3 = p / (3*u3) - u3 - a / 3.0

    if real:
        return [x.real
                for x in [root1, root2, root3]
                if abs(x.imag) < 1e-10]
    else:
        return [root1, root2, root3]


def _solveCubic_test(n=100):

    def test(a, b, c):
        xs = solveCubic(a, b, c)

        for x in xs:
            y = x**3 + a*x*x + b*x + c
            assert abs(y) < 1e-4, y

    test(0, 0, 0)
    test(0, 1, 1)
    test(0, 0, 1)

    for i in range(n):

        a = random.normalvariate(10, 5)
        b = random.normalvariate(10, 5)
        c = random.normalvariate(10, 5)

        test(a, b, c)




#=============================================================================
# testing

if __name__ == "__main__":

    # iter_window
    vals = sorted([random.random() * 20 for x in range(600)])

    vals += sorted([40 + random.random() * 20 for x in range(600)])

    '''
    win = filter(lambda x: len(x) > 0,
                 list(iter_window_index(vals, 5)))

    p = util.plot(util.cget(win, 2))#, style="lines")
    p.enableOutput(False)
    p.plot(util.cget(win, 3)) #, style="lines")

    for i, y in enumerate(vals):
        p.plot([i, len(vals)], [y, y], style="lines")
    p.enableOutput(True)
    p.replot()
    '''

    def mean2(v):
        if len(v) == 0:
            return 0.0
        else:
            return mean(v)

    x, y = zip(* iter_window_step(vals, 5, 1, len))
    # plotting removed (no gnuplot)
    # util.plot(x, y)
