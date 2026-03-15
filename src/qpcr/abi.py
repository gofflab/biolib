#!/usr/bin/env python
'''
Utilities for parsing and analyzing ABI qPCR instrument output.

Provides functions for parsing raw ABI results and cycle data files,
computing PCR amplification efficiencies via a sliding-window linear
regression on log-transformed fluorescence values, performing the
delta-delta Ct (ddCt) relative-quantification calculation, and
summarizing/reporting the results.

Requirements:
    - numpy
    - rpy
    - R (obviously)
        - lattice package (for plotting)

results.txt input format example (tab-delimited):
Well    Sample    Detector    Task    Ct    Threshold
1    cDNA_1    GapDH    EndogenousControl    25.678    0.6443

cycleData.txt input format example (tab-delimited):
Well    Sample    Detector    1    2    3    ...    nCycles
1    cDNA_1    GapDH    0.11    0.12    0.12    ...    6.57

Usage:
python abi.py results.txt cycleData.txt endoControl reference outFile

#TODO: change outFile to outDir

@author: lgoff
'''
###########################
#Imports
###########################
import math
import subprocess
import sys

import numpy as np

#from seqtools.misc import pp
#from rpy import *

###########################
#Constants
###########################
windowSize = 4
dictKeys = ['well','sample','detector','task','Ct','threshold']

##########################
#Parsing
##########################

def parseData(fname):
    """Parse a simplified ABI results text file into a list of well dictionaries.

    Raw input is a tab-delimited matrix with columns:
    Well, SampleName, DetectorName, Task, Ct, Threshold.
    Intermediate headers and summary rows must be removed from the raw ABI
    output before calling this function; only one header row should remain.
    Wells with an ``Undetermined`` Ct value are silently skipped.

    Args:
        fname: Path to the tab-delimited results text file.

    Returns:
        A list of dicts, one per well, with keys ``well`` (int),
        ``sample``, ``detector``, ``task``, ``Ct`` (float), and
        ``threshold`` (float).
    """
    data = []
    handle = open(fname,'r')
    #Remove Header Row
    headerRow = next(handle)
    headerVals = headerRow.rstrip().split('\t')
    #Parse well information
    for line in handle:
        values = line.rstrip().split('\t')
        if values[4] == 'Undetermined':
            continue
        myWell = dict(zip(dictKeys,values))
        myWell['well'],myWell['Ct'],myWell['threshold'] = int(myWell['well']),float(myWell['Ct']),float(myWell['threshold'])
        data.append(myWell)
    return data

def getDetAndSamp(data):
    """Return ordered lists of unique detector and sample names found in the data.

    Preserves first-seen order for both detectors and samples.

    Args:
        data: List of well dicts as returned by ``parseData``, each containing
            ``detector`` and ``sample`` keys.

    Returns:
        A tuple ``(detectors, samples)`` where each element is a list of
        unique string names in the order they were first encountered.
    """
    detectors = []
    samples = []
    for well in data:
        if well['detector'] not in detectors:
            detectors.append(well['detector'])
        if well['sample'] not in samples:
            samples.append(well['sample'])
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

def parseCycleData(fname):
    """Parse a tab-delimited cycle fluorescence file into a list of well dicts.

    Raw input is a tab-delimited file with a header row. Columns are:
    Well, Sample, Detector, followed by one column per cycle number.

    Args:
        fname: Path to the tab-delimited cycle data text file.

    Returns:
        A list of dicts, one per well, with keys ``well`` (int),
        ``sample`` (str), ``detector`` (str), and ``values`` (numpy array
        of float fluorescence readings, one per cycle).
    """
    cycleData = []
    handle = open(fname,'r')
    headerRow = next(handle)
    headerVals = headerRow.rstrip().split('\t')
    cycles = headerVals[3:]
    cycles = list(map(int,cycles))
    ncycles = int(headerVals[-1])

    for line in handle:
        values = line.rstrip().split('\t')
        well = int(values.pop(0))
        sample = values.pop(0)
        detector = values.pop(0)
        values = np.array(list(map(float,values)))
        cycleData.append({'well':well,'sample':sample, 'detector':detector, 'values': values})

    return cycleData

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
    median Ct for each combination.

    Args:
        data: List of well dicts, each containing ``sample``, ``detector``,
            and ``Ct`` keys.

    Returns:
        A nested dict ``{sample: {detector: median_Ct}}`` where each value
        is the median Ct (float) computed from all replicate wells.
    """
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

def calculateEfficiencies(cycleData):
    """Compute PCR amplification efficiency and initial concentration (N0) for each well.

    For each well, log10-transforms the fluorescence values, then slides a
    window of size ``windowSize`` across all cycles and picks the window with
    the highest Pearson correlation between log-fluorescence and cycle number
    (i.e., the most linear exponential-phase segment). A linear regression on
    that best window gives the slope (from which efficiency = 10^slope) and
    intercept (from which N0 = 10^intercept).

    Adds the following keys to each well dict in-place:
        ``logVals``, ``bestIdx``, ``bestCorr``, ``bestSlice``,
        ``bestCycles``, ``bestSlope``, ``bestIntercept``,
        ``efficiency``, ``N0``.

    Args:
        cycleData: List of well dicts as returned by ``parseCycleData``,
            each containing at minimum a ``values`` numpy array.

    Returns:
        The same list of well dicts with the additional efficiency and N0
        keys populated.
    """
    res = []
    for well in cycleData:
        well['logVals'] = getLogVals(well['values'])
        #Test all windows for linearity
        corrs = np.zeros(len(well['logVals'])-windowSize)
        for i in range(0,len(well['logVals'])-windowSize):
            logSlice = well['logVals'][i:i+windowSize]
            corrs[i]=corr(logSlice,np.array(range(1,windowSize+1)))
        #Append best Correlation Index to well
        well['bestIdx'] = np.argmax(corrs)

        #Do math on best window
        well['bestCorr'] = corrs[well['bestIdx']]
        well['bestSlice'] = np.array(well['logVals'][well['bestIdx']:well['bestIdx']+windowSize])
        well['bestCycles'] = np.array(range(well['bestIdx']+1,well['bestIdx']+1+windowSize))

        well['bestSlope'] = slope(well['bestCycles'],well['bestSlice'])
        well['bestIntercept'] = intercept(well['bestCycles'],well['bestSlice'])
        well['efficiency'] = 10**well['bestSlope']
        well['N0'] = 10**well['bestIntercept']
        res.append(well)
    return res

def summarizeEfficiencies(cycleData):
    """Compute mean and standard deviation of PCR efficiency for each detector.

    Groups per-well efficiency values by detector name and computes summary
    statistics.

    Args:
        cycleData: List of well dicts, each containing ``detector`` and
            ``efficiency`` keys (as produced by ``calculateEfficiencies``).

    Returns:
        A dict ``{detector: {'meanEff': float, 'sdevEff': float}}`` giving
        the mean and standard deviation of efficiency across all wells for
        each detector.
    """
    tmp = {}
    #Aggregate efficiencies by detector
    for i in cycleData:
        tmp.setdefault(i['detector'],[])
        tmp[i['detector']].append(i['efficiency'])
    #Get mean and sdev of efficiences by detector
    eff = {}
    for k in tmp.keys():
        eff[k] = {'meanEff':mean(tmp[k]),'sdevEff':sdev(tmp[k])}
    return eff

def mergeDataAndCycleData(data,cycleData,idx):
    """Copy efficiency and N0 values from cycleData into the matching well dicts in data.

    Uses the provided well-number index to look up each cycle-data well in
    the data list and transfers the ``N0`` and ``efficiency`` values. Wells
    present in cycleData but absent from data (e.g., wells skipped due to
    undetermined Ct) are silently ignored.

    Args:
        data: List of well dicts as returned by ``parseData``.
        cycleData: List of well dicts as returned by ``calculateEfficiencies``,
            each containing ``well``, ``N0``, and ``efficiency`` keys.
        idx: List of integer well numbers parallel to ``data``, as returned
            by ``wellIndex``.

    Returns:
        The ``data`` list with ``N0`` and ``efficiency`` keys added to each
        matched well dict.
    """
    for c in cycleData:
        try:
            dataloc = idx.index(c['well'])
            data[dataloc]['N0'] = c['N0']
            data[dataloc]['efficiency'] = c['efficiency']
        except ValueError:
            continue
    return data

#TODO: Make summarizer for N0 elements by sample and detector

def getLogVals(myArray):
    """Return the base-10 logarithm of each element in a numpy array.

    Args:
        myArray: A numpy array of positive numeric values.

    Returns:
        A numpy array of the same shape containing log10 of each input value.
    """
    return np.log10(myArray)

###############################
#ddCt math
###############################
def ddCt(data,medianCts,endoControl,reference):
    """Compute delta-Ct and delta-delta-Ct values for each well.

    For each well, dCt is calculated as:
        dCt = Ct - median_Ct(sample, endoControl)

    ddCt is then calculated as:
        ddCt = dCt - median_dCt(reference, detector)

    Wells where the endogenous control Ct is unavailable receive ``"N/A"``
    for dCt, and wells where the reference dCt is unavailable receive
    ``"N/A"`` for ddCt.

    Args:
        data: List of well dicts, each containing ``sample``, ``detector``,
            and ``Ct`` keys.
        medianCts: Nested dict ``{sample: {detector: median_Ct}}`` as returned
            by ``aggregateReplicateCts``.
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
        data[i]['dCt'] = data[i]['Ct'] - medianCts[data[i]['sample']][endoControl]
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
        except KeyError:
            data[i]['ddCt'] = "N/A"
            #print "%d\t%s" % (data[i]['well'],data[i]['ddCt'])
    return data

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
    """Compute the median of a list of numbers.

    Sorts the list in-place before computing.

    Args:
        vals: A list of numeric values.

    Returns:
        The median value as a float. For even-length lists, returns the
        average of the two middle values.
    """
    lenvals = len(vals)
    vals.sort()

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
    (i.e., one or both lists have zero variance), which is used as a
    sentinel for a perfect linear relationship in the sliding-window search.

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
    Wells with ``"N/A"`` RQ are excluded from RQ and dCt summaries but N0
    is always summarized (N0 values are assumed to be present).

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

def doPlotting(plotScript = "plotting.q"):
    """Execute an external R plotting script as a subprocess.

    Args:
        plotScript: Path to the R script to execute. Defaults to
            ``"plotting.q"``.

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
    """Run the full ABI qPCR analysis pipeline interactively.

    Parses results and cycle-data files, computes efficiencies,
    interactively asks the user to select an endogenous control and reference
    sample, performs ddCt/RQ calculations, and writes ``output.txt`` before
    running the external plotting script.

    Args:
        mainFile: Path to the tab-delimited ABI results file.
        cycleFile: Path to the tab-delimited cycle fluorescence file.
    """
    #Parse mainFile
    print("Parsing Results File...")
    data = parseData(mainFile)
    medianCts = aggregateReplicateCts(data) #Returns a dictionary of dictionaries by sample and then detector
    myIdx = wellIndex(data)

    #Efficiency Calculation from cycleFile
    print("Parsing CycleData File...")
    cycleData = parseCycleData(cycleFile)
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
