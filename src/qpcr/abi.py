#!/usr/bin/env python
'''
Created on Feb 22, 2010

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
import sys
import math
import numpy as np
import commands
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
    """Raw input for this file is a matrix of well x (Well,SampleName,DetectorName,Task,Ct,Threshold).  You must also delete the intermediate headers and summary rows from raw output of ABI.
    Be sure to remove the header section (except one header row).
    """
    data = []
    handle = open(fname,'r')
    #Remove Header Row
    headerRow = handle.next()
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
    detectors = []
    samples = []
    for well in data:
        if not well['detector'] in detectors:
            detectors.append(well['detector'])
        if not well['sample'] in samples:
            samples.append(well['sample'])
    return detectors,samples
    
def wellIndex(data):
    index = []
    for i in range(len(data)):
        index.append(data[i]['well'])
    return index

def parseCycleData(fname):
    """Raw input is tab-delimited text file with matrix of WellsxCycle values.  Header row is included.
    """
    cycleData = []
    handle = open(fname,'r')
    headerRow = handle.next()
    headerVals = headerRow.rstrip().split('\t')
    cycles = headerVals[3:]
    cycles = map(int,cycles)
    ncycles = int(headerVals[-1])
    
    for line in handle:
        values = line.rstrip().split('\t')
        well = int(values.pop(0))
        sample = values.pop(0)
        detector = values.pop(0)
        values = np.array(map(float,values))
        cycleData.append({'well':well,'sample':sample, 'detector':detector, 'values': values})
    
    return cycleData

######################
#Get User Input
######################
def getEndoControl(detectors):
    myString = "Please choose an endogenous control:\n"
    for i in range(0,len(detectors)):
        myString = myString+"\t(%d):\t%s\n" % (i,detectors[i])
    myString = myString + "Choose %s-%s:" % (0,len(detectors))
    choice = int(raw_input(myString))
    return detectors[choice]

def getReference(samples):
    myString = "Please choose a reference sample:\n"
    for i in range(0,len(samples)):
        myString = myString + "\t(%d):\t%s\n" % (i,samples[i])
    myString = myString + "Choose %s-%s:" % (0,len(samples))
    choice = int(raw_input(myString))
    return samples[choice]

#####################################
#Aggregate Replicates
#####################################

def aggregateReplicateCts(data):
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
    """Takes a list of dictionaries of cycle information by well and returns those same dictionaries with 
    additional keys for efficiency and concentration (N0) values."""
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
    """Takes an index of data (by well) and the cycleData to add the efficiency and N0 from cycleData to the 
    data dictionaries"""
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
    return np.log10(myArray)

###############################
#ddCt math
###############################
def ddCt(data,medianCts,endoControl,reference):
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
    """Computes the mean of a list of numbers"""
    n = 0
    s = 0.0
    for i in vals:
        s += i
        n += 1
    return s / float(n)

def median(vals):
    """Computes the median of a list of numbers"""
    lenvals = len(vals)
    vals.sort()
    
    if lenvals % 2 == 0:
        return (vals[lenvals / 2] + vals[lenvals / 2 - 1]) / 2.0
    else:
        return vals[lenvals / 2]

def variance(vals):
    """Variance"""
    u = mean(vals)
    return sum((x - u)**2 for x in vals) / float(len(vals)-1)

def sdev(vals):
    """Standard deviation"""
    if len(vals) <=1: return 0.0
    return math.sqrt(variance(vals))

def covariance(lst1, lst2):
    """Covariance"""
    m1 = mean(lst1)
    m2 = mean(lst2)
    tot = 0.0
    for i in xrange(len(lst1)):
        tot += (lst1[i] - m1) * (lst2[i] - m2)
    return tot / (len(lst1)-1)

def corr(lst1, lst2):
    """Pearson's Correlation"""
    num = covariance(lst1, lst2)
    denom = float(sdev(lst1) * sdev(lst2))
    if denom != 0:
        return num / denom
    else:
        return 1e1000

def slope(xarray,yarray):
    """Uses numpy, in fact assumes that the list arguments are numpy arrays."""
    n = float(len(xarray))
    m = (n*sum(xarray*yarray)-sum(xarray)*sum(yarray))/(n*sum(xarray**2)-(sum(xarray))**2)
    return m

def intercept(xarray,yarray):
    """Uses numpy, in fact assumes that the list arguments are numpy arrays."""
    m = slope(xarray,yarray)
    n = float(len(xarray))
    b = (sum(yarray)-m*(sum(xarray)))/n
    return b

###############################
#Reporting
###############################

def flagBadDetectors():
    pass

def aggregateResults(data):
    try:
        data[0]['RQ']
    except KeyError:
        print "Tried to aggregate RQs before they exist"
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
    #Open out Handle
    outHandle = open(outFile,'w')
    #Print header row
    print "Sample\tDetector\tmeanEff\tmeanRQ\tsdevRQ\tmedianRQ\tmeandCt\tmediandCt\tsdevdCt\tquant\tci.l\tci.u"
    print >>outHandle, "Sample\tDetector\tmeanEff\tmeanRQ\tsdevRQ\tmedianRQ\tmeandCt\tmediandCt\tsdevdCt\tquant\tci.l\tci.u"
    for sample,v in RQsummary.iteritems():
        for detector,v2 in v.iteritems():
            #print "%s\t%s\t%.2f\t%.2f\t%.2f" % (sample,detector,v2['meanRQ'],v2['medianRQ'],v2['sdevRQ'])
            print "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (sample,detector,effs[detector]['meanEff'],v2['meanRQ'],v2['sdevRQ'],v2['medianRQ'],v2['meandCt'],v2['mediandCt'],v2['sdevdCt'],effs[detector]['meanEff']**-v2['mediandCt'],effs[detector]['meanEff']**-(v2['mediandCt']+v2['sdevdCt']),effs[detector]['meanEff']**-(v2['mediandCt']-v2['sdevdCt']))
            print >>outHandle, "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (sample,detector,effs[detector]['meanEff'],v2['meanRQ'],v2['sdevRQ'],v2['medianRQ'],v2['meandCt'],v2['mediandCt'],v2['sdevdCt'],effs[detector]['meanEff']**-v2['mediandCt'],effs[detector]['meanEff']**-(v2['mediandCt']+v2['sdevdCt']),effs[detector]['meanEff']**-(v2['mediandCt']-v2['sdevdCt']))
    outHandle.close()

#######################
#Plotting
#######################
#TODO:Create R Function to plot output from printDataFramRQs()

def plotRQs(results):
    pass

def plotEdCt(results):
    pass

def doPlotting(plotScript = "plotting.q"):
    return commands.getstatusoutput(plotScript)
     

def makeDvsS(results,detectors,samples,value = "mediandCt"):
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
    #Parse mainFile
    print "Parsing Results File..."
    data = parseData(mainFile)
    medianCts = aggregateReplicateCts(data) #Returns a dictionary of dictionaries by sample and then detector
    myIdx = wellIndex(data)
    
    #Efficiency Calculation from cycleFile
    print "Parsing CycleData File..."
    cycleData = parseCycleData(cycleFile)
    cycleData = calculateEfficiencies(cycleData)
    effs = summarizeEfficiencies(cycleData)
    
    detectors,samples = getDetAndSamp(data)
    print "Found %d detectors (primers)..." % len(detectors)
    endoControl = getEndoControl(detectors)
    print "Found %d samples..." % len(samples)
    reference = getReference(samples)
    
    #Begin E^-ddCt Calculation
    data = ddCt(data,medianCts,endoControl,reference)
    data = RQ(data,effs)
    
    #Add effs and N0 from cycleData to well data
    data = mergeDataAndCycleData(data,cycleData,myIdx)
    
    #detectors,samples = getDetAndSamp(data)
    
    results = aggregateResults(data)
    printDataFrameRQs(results,effs,'output.txt')
    print "Output in 'output.txt'..."
    print "Plotting..."
    status = doPlotting()
    
    return
    
def test():
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
