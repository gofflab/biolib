#!/usr/bin/env python
'''
Core qPCR analysis module using four-parameter logistic modelling and iterative
nonlinear regression for efficiency estimation.

Provides the ``Well`` class for per-well data storage and curve fitting, along
with standalone functions for parsing raw ABI instrument output, performing
delta-delta Ct (ddCt) relative quantification, and reporting results.

This module extends the functionality in ``abi.py`` with a more rigorous
curve-fitting approach based on the four-parameter logistic (4PL) model
described in Zhao et al.

Requirements:
    - numpy
    - scipy

results.txt input format example (tab-delimited):
Well    Sample      Detector      Task      Ct    Threshold
1     cDNA_1       GapDH    EndogenousControl     25.678       0.6443

cycleData.txt input format example (tab-delimited):
Well    Sample      Detector      1       2    3     ...    nCycles
1     cDNA_1       GapDH    0.11    0.12    0.12    ...       6.57

Usage:
python abi.py results.txt cycleData.txt endoControl reference outFile

#TODO: change outFile to outDir

@author: lgoff
'''
###########################
#Imports
###########################
import itertools
import math
import subprocess
import sys

import numpy as np
from scipy import optimize

from . import util

#from seqtools.misc import pp
#from rpy import *

###########################
#Constants
###########################
windowSize = 4
dictKeys = ['well','sample','detector','task','Ct','threshold']

##########################
#Classes
##########################
class Well:
    """Represents a single PCR well with its raw data and fitted curve parameters.

    Stores metadata (sample name, detector, task, etc.), raw fluorescence
    readings keyed by cycle, and all intermediate and final results from
    four-parameter logistic curve fitting and crossing-point estimation.

    Attributes:
        wellNum: Integer well number (defaults to -1 until populated).
        sample: Sample name string.
        detector: Detector (primer/probe) name string.
        reporter: Reporter dye name string.
        task: Task type string (e.g., ``"EndogenousControl"``).
        Ct: Threshold cycle value (float).
        quantity: Quantity value from ABI output (float).
        eff: Amplification efficiency (float).
        threshold: Fluorescence threshold (float).
        cycles: List of cycle labels from the cycle data file.
        fluorData: Numpy array of fluorescence readings per cycle.
        flags: Dict of quality-flag name/value pairs parsed from the ABI file.
        RNoise: Standard error of the baseline fluorescence parameter (y0)
            from the fitted 4PL model; None until ``fitPCRCurve`` is called.
    """

    def __init__(self,line):
        """Initialise a Well with default empty values.

        Args:
            line: The raw text line from the ABI file used to create this
                well (stored for reference but not parsed here; parsing is
                done by ``parseRawABI``).
        """
        self.wellNum = -1
        self.sample = ''
        self.detector = ''
        self.reporter = ''
        self.task = ''
        self.Ct = 0.0
        self.quantity = 0.0
        self.eff = 0.0
        self.threshold = 0.0
        self.cycles = []
        self.fluorData = []
        self.flags = {}
        self.RNoise = None

    def estimateParams(self):
        """Generate initial parameter guesses for the four-parameter logistic model.

        Estimates starting values for the curve-fitting routine based on
        simple statistics of the raw fluorescence data:

        - ``y0``: mean of the first five cycles (baseline fluorescence).
        - ``x0``: cycle nearest the midpoint fluorescence (inflection point).
        - ``a``: dynamic range (max minus min fluorescence).
        - ``b``: set to 0 (the optimiser handles this parameter well without
          a manual initial estimate).

        Populates the instance attributes ``y0``, ``x0``, ``a``, and ``b``
        in-place.
        """
        self.y0 = np.mean(self.fluorData[:5]) # Initial guess as to baseline fluorescence (mean of first five cycles)
        self.x0 = self.cycles[np.argmin(abs(self.fluorData-np.mean(self.fluorData)))] # Initial guess as to inflection point at middle of curve
        self.a = (np.max(self.fluorData)-np.min(self.fluorData))# Initial guess as to y value at inflection
        self.b = 0 # Don't think I need to estimate this parameter, model seems to do a good job of fitting this one.

    def fitPCRCurve(self):
        """Fit the four-parameter logistic (4PL) model to the fluorescence data.

        Calls ``scipy.optimize.curve_fit`` with ``qpcrFit`` as the model
        function and up to 5000 function evaluations. After fitting,
        updates the instance attributes:

        - ``a``, ``b``, ``x0``, ``y0``: fitted model parameters.
        - ``pCov``: covariance matrix of the fitted parameters.
        - ``fitData``: list of model-predicted fluorescence values at each
          cycle.
        - ``paramSE``: dict mapping parameter names (``'a'``, ``'b'``,
          ``'x0'``, ``'y0'``) to their standard errors (sqrt of the
          diagonal of ``pCov``).
        - ``RNoise``: standard error of the ``y0`` parameter, used as an
          estimate of baseline noise.
        """
        #Fit qpcr Model
        newParams,self.pCov = optimize.curve_fit(qpcrFit,xdata=self.cycles,ydata=self.fluorData,maxfev=5000)
        #Update params
        self.a,self.b,self.x0,self.y0 = newParams
        #Generate fit data
        self.fitData = [qpcrFit(x,self.a,self.b,self.x0,self.y0) for x in self.cycles]
        #Find standard error of regression parameters as sqrt of variance from pCov
        self.paramSE = {}
        paramOrder = ['a','b','x0','y0']
        for i in range(4):
            self.paramSE[paramOrder[i]]=np.sqrt(self.pCov[i][i])
        #Get RNoise
        self.RNoise = self.paramSE['y0']
        return

    def CP_FDM(self):
        """Compute the crossing-point by the First Derivative Maximum (FDM) method.

        Calculates the cycle number at which the first derivative of the
        fitted 4PL curve is maximised, stored in ``self.FDM``.

        Returns:
            The FDM crossing-point cycle number as a float.
        """
        self.FDM = (self.x0*nthRoot(((self.b-1)/(self.b+1)),self.b))
        return self.FDM

    def CP_SDM(self):
        """Compute the crossing-point by the Second Derivative Maximum (SDM) method.

        Calculates the cycle number at which the second derivative of the
        fitted 4PL curve is maximised, stored in ``self.SDM``.

        Returns:
            The SDM crossing-point cycle number as a float.
        """
        self.SDM = self.x0*nthRoot((np.sqrt((3*self.b**2)*(self.b**2-1))-(2*(1-self.b**2)))/((self.b**2)+(3*self.b)+2),self.b)
        return self.SDM

    def CP_SPE(self):
        """Compute the crossing-point by the Signal-to-Noise method (SPE).

        Calculates the cycle at which the fluorescence signal exceeds the
        baseline noise by a factor of ``a / RNoise``, stored in ``self.SPE``.
        Requires that ``fitPCRCurve`` has been called so that ``RNoise`` is
        available.

        Returns:
            The SPE crossing-point cycle number as a float.
        """
        self.SPE = (self.x0*nthRoot(((self.a-self.RNoise)/self.RNoise),self.b))
        return self.SPE

    def iterativeNLR(self):
        """Perform iterative nonlinear regression over the exponential phase window.

        Uses the SPE and SDM crossing-point estimates to define the lower and
        upper cycle boundaries of the exponential phase. Enumerates all
        sub-windows of size >= ``windowSize`` within that range using
        combinatorics and stores the window indices in ``self.winIdx``.

        Requires that ``CP_SPE`` and ``CP_SDM`` have been called first to
        populate ``self.SPE`` and ``self.SDM``.
        """
        self.lowerCycleNum = int(self.SPE)
        self.upperCycleNum = int(self.SDM)
        self.regPoints = self.upperCycleNum-self.lowerCycleNum+1
        #Get windows
        winIdx = []
        for i in range(windowSize,self.regPoints+1):
            combs = itertools.combinations(range(self.lowerCycleNum,self.upperCycleNum+1),i)
            for c in combs:
                winIdx.append(c)



##########################
#Parsing
##########################
def parseRawABI(fname):
    """Parse a raw ABI results file into a dict of Well objects keyed by well number.

    Replaces the simpler ``parseData`` function. Handles the multi-section ABI
    export format: skips the first line, collects key/value header metadata,
    then reads data rows until EOF. Rows with an ``"Undetermined"`` Ct value
    are skipped. Quality-flag columns (indices 17 onwards) are stored in each
    ``Well.flags`` dict.

    Args:
        fname: Path to the raw ABI tab-delimited results export file.

    Returns:
        A dict ``{well_number (int): Well}`` for every well with a valid
        numeric Ct value.
    """
    dictKeys = ['well','sample','detector','reporter','task','Ct','quantity','Qty Mean','Qty StdDev','Ct Median','Ct Mean','Ct StdDev','Baseline Type','Baseline Start','Baseline Stop','Threshold Type','threshold','FOS','HMD','LME','EW','BPR','NAW','HNS','HRN','EAF','BAF','TAF','CAF']
    handle = open(fname,'r')
    header = {}
    res = {}
    handle.readline()#Skip first line

    #Collect header information
    while True:
        line = handle.readline()
        if line.startswith("Well"):
            break
        vals = line.rstrip("\r\n").split("\t")
        if len(vals)==2:
            header[vals[0]]=vals[1]

    while True:
        if line.startswith("Well"):
            #print line
            #myKeys = line.rstrip("\r\n").split("\t")[:-1]
            line = handle.readline()
            continue
        vals = line.rstrip("\r\n").split("\t")[:-1]
        if len(vals)!=len(dictKeys):
            line=handle.readline()
            if not line: break
            continue
        vals[0] = int(vals[0])
        #Ignore anything that has an undetermined Ct value
        if vals[5] == "Undetermined":
            line = handle.readline()
            continue

        for i in [5,6,7]:
            try:
                vals[i] = float(vals[i])
            except ValueError:
                pass
        try:
            vals[16] = float(vals[16])
        except ValueError:
            pass
        try:
            tmp = dict(zip(dictKeys,vals))
            myWell = Well(line)
            myWell.wellNum,myWell.sample,myWell.detector,myWell.reporter,myWell.task,myWell.threshold,myWell.flags = tmp['well'],tmp['sample'],tmp['detector'],tmp['reporter'],tmp['task'],tmp['threshold'],dict(zip(dictKeys[17:],vals[17:]))
            res[myWell.wellNum] = myWell
        except ValueError:
            pass
        line=handle.readline()
        if not line: break
    return res

    assert False, "Should not reach this line..."

def parseRawCycle(fname,wellData):
    """Parse a raw ABI cycle fluorescence file and populate the matching Well objects.

    Replaces the simpler ``parseCycleData`` function. Reads fluorescence
    readings up to (but not including) the ``"Delta Rn"`` column and writes
    ``cycles`` and ``fluorData`` directly onto the corresponding ``Well``
    objects in ``wellData``.

    Args:
        fname: Path to the raw ABI cycle data tab-delimited export file.
        wellData: Dict ``{well_number: Well}`` as returned by
            ``parseRawABI``. Modified in-place.
    """
    handle = open(fname,'r')
    handle.readline()#Remove first line
    headerRow = handle.readline()
    headerVals = headerRow.rstrip().split("\t")
    myLim = headerVals.index("Delta Rn")
    headerVals = headerVals[:myLim]
    for line in handle:
        vals = line.rstrip().split("\t")[:myLim]
        well = int(vals.pop(0))
        detector = vals.pop(0)
        vals = np.array(list(map(float,vals[1:])))
        wellData[well].cycles,wellData[well].fluorData = headerVals,vals
    return

def getDetAndSamp(wellData):
    """Return lists of unique detector and sample names from a collection of Well objects.

    Uses ``util.uniqify`` to deduplicate; result order is not guaranteed to
    be preserved (depends on dict key ordering).

    Args:
        wellData: An iterable of ``Well`` objects (e.g., the values of the
            dict returned by ``parseRawABI``).

    Returns:
        A tuple ``(detectors, samples)`` where each element is a list of
        unique string names.
    """
    detectors = util.uniqify(detectors = [x.detector for x in wellData])
    samples = util.uniqify(samples = [x.sample for x in wellData])
    return detectors,samples

def wellIndex(data):
    """Build a list of well numbers in the same order as the data list.

    Args:
        data: List of well dicts, each containing a ``well`` key.

    Returns:
        A list of integer well numbers corresponding positionally to each
        entry in ``data``.
    """
    index = []
    for i in range(len(data)):
        index.append(data[i]['well'])
    return index

######################
#Get User Input
######################
def getEndoControl(detectors):
    """Interactively prompt the user to select an endogenous control detector.

    Prints a numbered list of detector names and reads an integer choice from
    standard input.

    Args:
        detectors: List of detector name strings to present to the user.

    Returns:
        The detector name string chosen by the user.
    """
    myString = "Please choose an endogenous control:\n"
    for i in range(0,len(detectors)):
        myString = myString+"\t(%d):\t%s\n" % (i,detectors[i])
    myString = myString + "Choose %s-%s:" % (0,len(detectors))
    choice = int(input(myString))
    return detectors[choice]

def getReference(samples):
    """Interactively prompt the user to select a reference sample.

    Prints a numbered list of sample names and reads an integer choice from
    standard input.

    Args:
        samples: List of sample name strings to present to the user.

    Returns:
        The sample name string chosen by the user.
    """
    myString = "Please choose a reference sample:\n"
    for i in range(0,len(samples)):
        myString = myString + "\t(%d):\t%s\n" % (i,samples[i])
    myString = myString + "Choose %s-%s:" % (0,len(samples))
    choice = int(input(myString))
    return samples[choice]

#####################################
#Aggregate Replicates
#####################################

def aggregateReplicateCts(data):
    """Aggregate replicate Ct values per sample/detector pair using the median.

    Groups raw per-well Ct values by (sample, detector) and computes the
    median Ct for each combination. ``"N/A"`` values (from undetermined wells
    that slipped through) are silently dropped by the ``median`` helper.

    Args:
        data: List of well dicts, each containing ``sample``, ``detector``,
            and ``Ct`` keys.

    Returns:
        A nested dict ``{sample: {detector: median_Ct}}`` where each value
        is the median Ct (float or ``"N/A"`` if all replicates are missing).
    """
    #This will have to change...
    #TODO: make this aggregate either Ct values or N0 values?
    tmp = {}
    for d in data:
        #print d
        tmp.setdefault(d['sample'],{})
        #print tmp
        tmp[d['sample']].setdefault(d['detector'],[])
        tmp[d['sample']][d['detector']].append(d['Ct'])
    medians = {}
    for k1 in tmp.keys():
        medians.setdefault(k1,{})
        for k2 in tmp[k1].keys():
            #print tmp[k1][k2]
            medians[k1][k2] = median(tmp[k1][k2])
    return medians

#####################################
#Calculate Efficiencies by Detector
#####################################

def getLogVals(myArray):
    """Return the base-10 logarithm of each element in a numpy array.

    Args:
        myArray: A numpy array of positive numeric values.

    Returns:
        A numpy array of the same shape containing log10 of each input value.
    """
    return np.log10(myArray)

#########
# Four-parameter Logistic Model fitting
#########
def nthRoot(num,n):
    """Compute the nth root of a number.

    Args:
        num: The base value (numeric).
        n: The root degree (numeric, must not be zero).

    Returns:
        ``num ** (1.0 / n)`` as a float.
    """
    return num ** (1.0/n)

def qpcrFit(x,a,b,x0,y0):
    """Evaluate the four-parameter logistic (4PL) model for qPCR fluorescence data.

    Implements the model from Zhao et al.:
        f(x) = y0 + a / (1 + (x / x0)^b)

    Designed for use with ``scipy.optimize.curve_fit``.

    Args:
        x: Cycle number (scalar or array).
        a: Amplitude parameter (difference between upper and lower
           asymptotes).
        b: Slope/steepness parameter.
        x0: Inflection point (cycle at the midpoint of the curve).
        y0: Baseline fluorescence (lower asymptote).

    Returns:
        Predicted fluorescence value(s) at cycle ``x``.
    """
    return (y0+(a/(1+((x/x0)**b))))

def qpcrFitResiduals(x,y,a,b,x0,y0):
    """Compute residuals between observed fluorescence and the 4PL model.

    Calculates ``y - qpcrFit(x, a, b, x0, y0)``.

    Args:
        x: Cycle number(s) (scalar or array).
        y: Observed fluorescence value(s).
        a: Amplitude parameter.
        b: Slope/steepness parameter.
        x0: Inflection point (cycle at midpoint).
        y0: Baseline fluorescence (lower asymptote).

    Returns:
        Residual value(s) ``y - predicted``.
    """
    return y-qpcrFit(x,a,b,x0,y0)

def CP_FDM(p):
    """Compute the crossing-point using the First Derivative Maximum (FDM) method.

    Args:
        p: Sequence of four fitted 4PL parameters ``[a, b, x0, y0]``.

    Returns:
        The FDM crossing-point cycle number as a float.
    """
    return (p[2]*nthRoot(((p[1]-1)/(p[1]+1)),p[1]))

def CP_SDM(p):
    """Compute the crossing-point using the Second Derivative Maximum (SDM) method.

    Args:
        p: Sequence of four fitted 4PL parameters ``[a, b, x0, y0]``.

    Returns:
        The SDM crossing-point cycle number as a float.
    """
    return p[2]*nthRoot((np.sqrt((3*p[1]**2)*(p[1]**2-1))-(2*(1-p[1]**2)))/((p[1]**2)+(3*p[1])+2),p[1])

def CP_SPE(p,rNoise):
    """Compute the crossing-point using the Signal-to-Noise (SPE) method.

    Args:
        p: Sequence of four fitted 4PL parameters ``[a, b, x0, y0]``.
        rNoise: Baseline noise estimate (standard error of the ``y0``
            parameter, i.e., ``RNoise``).

    Returns:
        The SPE crossing-point cycle number as a float.
    """
    return (p[2]*nthRoot(((p[0]-rNoise)/rNoise),p[1]))

###############################
#Iterative Nonlinear Regression
###############################
def nlmFit(x,a,b,y0):
    """Evaluate the exponential nonlinear regression model for the exponential phase.

    Models the exponential amplification phase as:
        f(x) = y0 + a * (b ^ x)

    Used for iterative nonlinear regression (iNLR) on windows within the
    exponential phase. Parameters are ``[a, b, y0]``.

    Args:
        x: Cycle number (scalar or array).
        a: Amplitude scaling factor.
        b: Per-cycle amplification factor (related to efficiency: b ~ E).
        y0: Baseline offset.

    Returns:
        Predicted fluorescence value(s) at cycle ``x``.
    """
    return y0+(a*(b**x))

def nlmFitResiduals(x,y,a,b,y0):
    """Compute residuals between observed fluorescence and the exponential NLM model.

    Calculates ``y - nlmFit(x, a, b, y0)``.

    Args:
        x: Cycle number(s) (scalar or array).
        y: Observed fluorescence value(s).
        a: Amplitude scaling factor.
        b: Per-cycle amplification factor.
        y0: Baseline offset.

    Returns:
        Residual value(s) ``y - predicted``.
    """
    return y-nlmFit(x,a,b,y0)


###############################
#ddCt math
###############################
def ddCt(data,medianCts,endoControl,reference):
    """Compute delta-Ct and delta-delta-Ct values for each well.

    For each well, dCt is calculated as:
        dCt = Ct - median_Ct(sample, endoControl)

    If the endogenous control Ct is unavailable for a sample, dCt is set to
    ``"N/A"``. ddCt is then calculated as:
        ddCt = dCt - median_dCt(reference, detector)

    If the reference dCt is unavailable, ddCt is set to ``"N/A"``.

    Args:
        data: List of well dicts, each containing ``sample``, ``detector``,
            and ``Ct`` keys.
        medianCts: Nested dict ``{sample: {detector: median_Ct}}`` as
            returned by ``aggregateReplicateCts``.
        endoControl: Name of the endogenous control detector to use for
            normalization.
        reference: Name of the reference sample to use for ddCt calculation.

    Returns:
        The ``data`` list with ``dCt`` and ``ddCt`` keys added to each well
        dict (values are floats or ``"N/A"``).
    """
    tmp = {}
    #Calculate dCts
    for i in range(len(data)):
        print(medianCts[data[i]['sample']])
        try:
            data[i]['dCt'] = data[i]['Ct'] - medianCts[data[i]['sample']][endoControl]
        except KeyError:
            data[i]['dCt'] = "N/A"
        tmp.setdefault(data[i]['sample'],{})
        tmp[data[i]['sample']].setdefault(data[i]['detector'],[])
        tmp[data[i]['sample']][data[i]['detector']].append(data[i]['dCt'])
    #Calculate median dCts
    med = {}
    for k1 in tmp.keys():
        med.setdefault(k1,{})
        for k2 in tmp[k1].keys():
            #print tmp[k1][k2]
            med[k1][k2] = median(tmp[k1][k2])

    #Calculate ddCts
    for i in range(len(data)):
        try:
            data[i]['ddCt'] = data[i]['dCt'] - med[reference][data[i]['detector']]
            #print "%d\t%.2f" % (data[i]['well'],data[i]['ddCt'])
        except:
            data[i]['ddCt'] = "N/A"
            #print "%d\t%s" % (data[i]['well'],data[i]['ddCt'])
    return data

def JohnsMethod(data,medianCts,endoControl,reference):
    """Placeholder for an alternative relative quantification method.

    Not yet implemented.

    Args:
        data: List of well dicts.
        medianCts: Nested dict of median Ct values per sample/detector.
        endoControl: Name of the endogenous control detector.
        reference: Name of the reference sample.
    """
    pass

def RQ(data,effs):
    """Calculate relative quantification (RQ) values for each well.

    RQ is computed as:
        RQ = meanEfficiency ^ (-ddCt)

    Wells with a ``"N/A"`` ddCt or a missing efficiency entry receive
    ``"N/A"`` for RQ.

    Args:
        data: List of well dicts containing ``detector`` and ``ddCt`` keys,
            as returned by ``ddCt``.
        effs: Dict ``{detector: {'meanEff': float, ...}}`` as returned by
            ``summarizeEfficiencies``.

    Returns:
        The ``data`` list with an ``RQ`` key added to each well dict
        (float or ``"N/A"``).
    """
    res = []
    for d in data:
        try:
            d['RQ'] = effs[d['detector']]['meanEff']**(-d['ddCt'])
        except:
            d['RQ'] = "N/A"
        res.append(d)
        #print "%d\t%s" % (d['well'],d['RQ'])
    return res



###############################
#Statistics
###############################

def mean(vals):
    """Compute the arithmetic mean of a list of numbers.

    Args:
        vals: An iterable of numeric values.

    Returns:
        The arithmetic mean as a float.
    """
    n = 0
    s = 0.0
    for i in vals:
        s += i
        n += 1
    return s / float(n)

def median(vals):
    """Compute the median of a list, ignoring any ``"N/A"`` sentinel values.

    Filters out ``"N/A"`` entries before sorting. Sorts the remaining values
    in-place. Returns ``"N/A"`` if no numeric values remain after filtering.

    Args:
        vals: A list that may contain numeric values and/or the string
            ``"N/A"``.

    Returns:
        The median numeric value as a float, or the string ``"N/A"`` if all
        values are ``"N/A"``.
    """
    print(vals)
    vals = [i for i in vals if i != "N/A"]
    print(vals)
    lenvals = len(vals)
    vals.sort()
    if lenvals == 0:
        return "N/A"
    if lenvals % 2 == 0:
        return (vals[lenvals // 2] + vals[lenvals // 2 - 1]) / 2.0
    else:
        return vals[lenvals // 2]

def variance(vals):
    """Compute the sample variance of a list of numbers.

    Uses Bessel's correction (divides by N-1).

    Args:
        vals: A list of numeric values with at least two elements.

    Returns:
        The sample variance as a float.
    """
    u = mean(vals)
    return sum((x - u)**2 for x in vals) / float(len(vals)-1)

def sdev(vals):
    """Compute the sample standard deviation of a list of numbers.

    Returns 0.0 for lists with one or fewer elements.

    Args:
        vals: A list of numeric values.

    Returns:
        The sample standard deviation as a float.
    """
    if len(vals) <=1: return 0.0
    return math.sqrt(variance(vals))

def covariance(lst1, lst2):
    """Compute the sample covariance between two equal-length lists.

    Uses Bessel's correction (divides by N-1).

    Args:
        lst1: First list of numeric values.
        lst2: Second list of numeric values; must be the same length as
            ``lst1``.

    Returns:
        The sample covariance as a float.
    """
    m1 = mean(lst1)
    m2 = mean(lst2)
    tot = 0.0
    for i in range(len(lst1)):
        tot += (lst1[i] - m1) * (lst2[i] - m2)
    return tot / (len(lst1)-1)

def corr(lst1, lst2):
    """Compute the Pearson correlation coefficient between two lists.

    Returns a very large number (1e1000) when the denominator is zero
    (i.e., one or both lists have zero variance), used as a sentinel for
    a perfect linear relationship in the sliding-window search.

    Args:
        lst1: First list of numeric values.
        lst2: Second list of numeric values; must be the same length as
            ``lst1``.

    Returns:
        The Pearson correlation coefficient as a float, or 1e1000 when the
        standard deviation of either list is zero.
    """
    num = covariance(lst1, lst2)
    denom = float(sdev(lst1) * sdev(lst2))
    if denom != 0:
        return num / denom
    else:
        return 1e1000

def slope(xarray,yarray):
    """Compute the ordinary least-squares regression slope.

    Uses the standard closed-form formula. Requires numpy arrays because
    element-wise multiplication (``xarray * yarray``) and vectorized
    ``sum`` are used.

    Args:
        xarray: Numpy array of independent variable values.
        yarray: Numpy array of dependent variable values; must be the same
            length as ``xarray``.

    Returns:
        The regression slope as a float.
    """
    n = float(len(xarray))
    m = (n*sum(xarray*yarray)-sum(xarray)*sum(yarray))/(n*sum(xarray**2)-(sum(xarray))**2)
    return m

def intercept(xarray,yarray):
    """Compute the ordinary least-squares regression intercept.

    Uses the standard closed-form formula given the slope. Requires numpy
    arrays because vectorized ``sum`` is used.

    Args:
        xarray: Numpy array of independent variable values.
        yarray: Numpy array of dependent variable values; must be the same
            length as ``xarray``.

    Returns:
        The regression intercept (y-axis) as a float.
    """
    m = slope(xarray,yarray)
    n = float(len(xarray))
    b = (sum(yarray)-m*(sum(xarray)))/n
    return b

###############################
#Reporting
###############################

def flagBadDetectors():
    """Flag detectors with poor amplification characteristics.

    Not yet implemented.
    """
    pass

def aggregateResults(data):
    """Aggregate per-well RQ, N0, and dCt values into per-(sample, detector) summaries.

    Computes mean, median, and standard deviation of RQ, dCt, and N0
    for every (sample, detector) combination across all replicate wells.
    Wells with ``"N/A"`` RQ are excluded from RQ and dCt summaries; N0
    is always summarised.

    Args:
        data: List of well dicts containing ``sample``, ``detector``,
            ``RQ``, ``N0``, and ``dCt`` keys, as returned by ``RQ``.

    Returns:
        A nested dict ``{sample: {detector: stats_dict}}`` where
        ``stats_dict`` contains the keys: ``medianRQ``, ``meanRQ``,
        ``sdevRQ``, ``mediandCt``, ``meandCt``, ``sdevdCt``,
        ``medianN0``, ``meanN0``, ``sdevN0``. Unavailable values are
        represented as ``"N/A"``.

    Raises:
        KeyError: If ``RQ`` values have not yet been computed on the data
            (i.e., ``ddCt`` and ``RQ`` have not been called first).
    """
    try:
        data[0]['RQ']
    except KeyError:
        print("Tried to aggregate RQs before they exist")
        raise
    #Setup intermediate lists to aggregate later
    tmpRQ = {}
    tmpN0 = {}
    tmpdCt = {}

    for d in data:
        if d['RQ'] == "N/A": continue
        #print d
        tmpRQ.setdefault(d['sample'],{})
        tmpN0.setdefault(d['sample'],{})
        tmpdCt.setdefault(d['sample'],{})
        #print tmp
        tmpRQ[d['sample']].setdefault(d['detector'],[])
        tmpN0[d['sample']].setdefault(d['detector'],[])
        tmpdCt[d['sample']].setdefault(d['detector'],[])

        tmpRQ[d['sample']][d['detector']].append(d['RQ'])
        tmpN0[d['sample']][d['detector']].append(d['N0'])
        tmpdCt[d['sample']][d['detector']].append(d['dCt'])

    #Aggregate temporary lists
    res = {}
    for k1 in tmpRQ.keys():
        res.setdefault(k1,{})
        for k2 in tmpRQ[k1].keys():
            #print tmp[k1][k2]
            res[k1].setdefault(k2,{})
            #Summarize RQ values
            RQlist = tmpRQ[k1][k2]
            naCount = RQlist.count("N/A")
            if naCount == len(RQlist):
                res[k1][k2]['medianRQ'] = "N/A"
                res[k1][k2]['meanRQ'] = "N/A"
                res[k1][k2]['sdevRQ'] = "N/A"

                res[k1][k2]['mediandCt'] = "N/A"
                res[k1][k2]['meandCt'] = "N/A"
                res[k1][k2]['sdevdCt'] = "N/A"
            else:
                for i in range(0,naCount):
                    RQlist.remove("N/A")
                res[k1][k2]['medianRQ'] = median(RQlist)
                res[k1][k2]['meanRQ'] = mean(RQlist)
                res[k1][k2]['sdevRQ'] = sdev(RQlist)

                    #Summarize dCt values
                res[k1][k2]['mediandCt'] = median(tmpdCt[k1][k2])
                res[k1][k2]['meandCt'] = mean(tmpdCt[k1][k2])
                res[k1][k2]['sdevdCt'] = sdev(tmpdCt[k1][k2])

            #Summarize N0 values (Possibly delete this later)
            res[k1][k2]['medianN0'] = median(tmpN0[k1][k2])
            res[k1][k2]['meanN0'] = mean(tmpN0[k1][k2])
            res[k1][k2]['sdevN0'] = sdev(tmpN0[k1][k2])

    return res

def printDataFrameRQs(RQsummary,effs,outFile):
    """Write a tab-delimited summary of RQ results to a file and to stdout.

    Outputs one row per (sample, detector) combination with columns:
    Sample, Detector, meanEff, meanRQ, sdevRQ, medianRQ, meandCt,
    mediandCt, sdevdCt, quant, ci.l, ci.u.

    The ``quant`` column is efficiency^(-mediandCt); ``ci.l`` and ``ci.u``
    are efficiency^(-(mediandCt +/- sdevdCt)), providing approximate
    confidence intervals.

    Args:
        RQsummary: Nested dict as returned by ``aggregateResults``.
        effs: Dict ``{detector: {'meanEff': float, ...}}`` as returned by
            ``summarizeEfficiencies``.
        outFile: Path to the output file to write.
    """
    #Open out Handle
    outHandle = open(outFile,'w')
    #Print header row
    print("Sample\tDetector\tmeanEff\tmeanRQ\tsdevRQ\tmedianRQ\tmeandCt\tmediandCt\tsdevdCt\tquant\tci.l\tci.u")
    print("Sample\tDetector\tmeanEff\tmeanRQ\tsdevRQ\tmedianRQ\tmeandCt\tmediandCt\tsdevdCt\tquant\tci.l\tci.u", file=outHandle)
    for sample,v in RQsummary.items():
        for detector,v2 in v.items():
            #print "%s\t%s\t%.2f\t%.2f\t%.2f" % (sample,detector,v2['meanRQ'],v2['medianRQ'],v2['sdevRQ'])
            print("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (sample,detector,effs[detector]['meanEff'],v2['meanRQ'],v2['sdevRQ'],v2['medianRQ'],v2['meandCt'],v2['mediandCt'],v2['sdevdCt'],effs[detector]['meanEff']**-v2['mediandCt'],effs[detector]['meanEff']**-(v2['mediandCt']+v2['sdevdCt']),effs[detector]['meanEff']**-(v2['mediandCt']-v2['sdevdCt'])))
            print("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (sample,detector,effs[detector]['meanEff'],v2['meanRQ'],v2['sdevRQ'],v2['medianRQ'],v2['meandCt'],v2['mediandCt'],v2['sdevdCt'],effs[detector]['meanEff']**-v2['mediandCt'],effs[detector]['meanEff']**-(v2['mediandCt']+v2['sdevdCt']),effs[detector]['meanEff']**-(v2['mediandCt']-v2['sdevdCt'])), file=outHandle)
    outHandle.close()

#######################
#Plotting
#######################
#TODO:Create R Function to plot output from printDataFramRQs()

def plotRQs(results):
    """Plot relative quantification (RQ) values.

    Not yet implemented.

    Args:
        results: Aggregated results dict as returned by ``aggregateResults``.
    """
    pass

def plotEdCt(results):
    """Plot efficiency-corrected delta-Ct (EdCt) values.

    Not yet implemented.

    Args:
        results: Aggregated results dict as returned by ``aggregateResults``.
    """
    pass

def doPlotting(plotScript = "qPCRPlotting.q"):
    """Execute an external R plotting script as a subprocess.

    Args:
        plotScript: Path to the R script to execute. Defaults to
            ``"qPCRPlotting.q"``.

    Returns:
        A tuple ``(status, output)`` as returned by
        ``subprocess.getstatusoutput``.
    """
    return subprocess.getstatusoutput(plotScript)


def makeDvsS(results,detectors,samples,value = "mediandCt"):
    """Build a detector-by-sample matrix of a chosen summary statistic.

    Creates a 2-D numpy array indexed by detector (rows) and sample
    (columns). Missing (sample, detector) combinations are filled with
    ``nan``.

    Args:
        results: Nested dict ``{sample: {detector: stats_dict}}`` as
            returned by ``aggregateResults``.
        detectors: Ordered list of detector names defining the row order.
        samples: Ordered list of sample names defining the column order.
        value: Key within the innermost stats dict to extract. Defaults to
            ``"mediandCt"``.

    Returns:
        A numpy float array of shape ``(len(detectors), len(samples))``
        containing the requested statistic for each cell.
    """
    matrix = np.zeros((len(detectors),len(samples)),float)
    for d in range(0,len(detectors)):
        for s in range(0,len(samples)):
            try:
                matrix[d,s]=results[samples[s]][detectors[d]][value]
            except KeyError:
                matrix[d,s]='nan'
    return matrix

##############################
#Main
##############################

def main(mainFile,cycleFile):
    """Run the full qPCR analysis pipeline interactively.

    Parses results and cycle-data files using the raw ABI format parsers,
    computes efficiencies, interactively asks the user to select an endogenous
    control and reference sample, performs ddCt/RQ calculations, writes
    ``output.txt``, and runs the external plotting script.

    Args:
        mainFile: Path to the raw ABI tab-delimited results export file.
        cycleFile: Path to the raw ABI cycle fluorescence export file.
    """
    #Parse mainFile
    print("Parsing Results File...")
    data = parseRawABI(mainFile)
    medianCts = aggregateReplicateCts(data) #Returns a dictionary of dictionaries by sample and then detector
    myIdx = wellIndex(data)

    #Efficiency Calculation from cycleFile
    print("Parsing CycleData File...")
    cycleData = parseRawCycle(cycleFile)
    cycleData = calculateEfficiencies(cycleData)
    effs = summarizeEfficiencies(cycleData)

    detectors,samples = getDetAndSamp(data)
    print("Found %d detectors (primers)..." % len(detectors))
    endoControl = getEndoControl(detectors)
    print("Found %d samples..." % len(samples))
    reference = getReference(samples)

    #Begin E^-ddCt Calculation
    data = ddCt(data,medianCts,endoControl,reference)
    data = RQ(data,effs)

    #Add effs and N0 from cycleData to well data
    data = mergeDataAndCycleData(data,cycleData,myIdx)

    #detectors,samples = getDetAndSamp(data)

    results = aggregateResults(data)
    printDataFrameRQs(results,effs,'output.txt')
    print("Output in 'output.txt'...")
    print("Plotting...")
    status = doPlotting()

    return

def test():
    """Run a manual integration test using hard-coded HeLa RIP data files.

    Parses ``'RIP HeLa clipped.txt'`` and ``'new_RIP_HeLa.txt'``, runs the
    full ddCt/RQ pipeline with hard-coded endogenous control (``'hGAPDH'``)
    and reference sample (``'IgG RIP'``), writes ``output.txt``, and
    returns a detector-by-sample matrix of mediandCt values.

    Returns:
        A numpy float array of shape ``(n_detectors, n_samples)`` containing
        the mediandCt for each (detector, sample) combination.
    """
    cycleData = parseCycleData('RIP HeLa clipped.txt')
    cycleData = calculateEfficiencies(cycleData)
    effs = summarizeEfficiencies(cycleData)
    #pp(effs)
    data = parseData('new_RIP_HeLa.txt')
    myIdx = wellIndex(data)
    summary = aggregateReplicateCts(data)
    #pp(summary)
    data = ddCt(data,summary,'hGAPDH','IgG RIP')
    #pp(data)
    data = RQ(data,effs)
    data = mergeDataAndCycleData(data,cycleData,myIdx)
    #pp(data)

    #Get Unique detectors and Sample Names to aid in plotting
    detectors,samples = getDetAndSamp(data)

    results = aggregateResults(data)
    #pp(results)
    printDataFrameRQs(results,effs,'output.txt')
    myMat = makeDvsS(results,detectors,samples)

    return myMat

if __name__ == '__main__':
    mainFile = sys.argv[1]
    cycleFile = sys.argv[2]
    main(mainFile,cycleFile)
