#!/usr/bin/env python
'''
Created on Feb 22, 2010

Requirements:
    - numpy
    - rpy
    - R (obviously)
        - lattice package (for plotting)

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
import sys
import math
import numpy as np
from scipy import optimize
import commands
import util
import itertools
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
    def __init__(self,line):
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
        self.y0 = np.mean(self.fluorData[:5]) # Initial guess as to baseline fluorescence (mean of first five cycles)
        self.x0 = self.cycles[np.argmin(abs(self.fluorData-np.mean(self.fluorData)))] # Initial guess as to inflection point at middle of curve
        self.a = (np.max(self.fluorData)-np.min(self.fluorData))# Initial guess as to y value at inflection 
        self.b = 0 # Don't think I need to estimate this parameter, model seems to do a good job of fitting this one.
    
    def fitPCRCurve(self):
        #Fit qpcr Model
        newParams,self.pCov = optimize.curvefit(qpcrFit,xdata=self.cycles,ydata=self.fluorData,maxfev=5000)
        #Update params
        self.a,self.b,self.x0,self.y0 = newParams
        #Generate fit data
        self.fitData = [qpcrFit(x,self.a,self.b,self.x0,self.y0) for x in self.cycles]
        #Find standard error of regression parameters as sqrt of variance from pCov
        self.paramSE = {}
        paramOrder = ['a','b','x0','y0']
        for i in xrange(4):
            self.paramSE[paramOrder[i]]=np.sqrt(self.pCov[i][i])
        #Get RNoise
        self.RNoise = self.paramSE['y0']
        return
    
    def CP_FDM(self):
        self.FDM = (self.x0*nthRoot(((self.b-1)/(self.b+1)),self.b))
        return self.FDM
    
    def CP_SDM(self):
        self.SDM = self.x0*nthRoot((np.sqrt((3*self.b**2)*(self.b**2-1))-(2*(1-self.b**2)))/((self.b**2)+(3*self.b)+2),self.b)
        return self.SDM
    
    def CP_SPE(self):
        self.SPE = (self.x0*nthRoot(((self.a-self.RNoise)/self.RNoise),self.b))
        return self.SPE
    
    def iterativeNLR(self):
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
    """This replaces parseData"""
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
            myWell = Well()
            myWell.wellNum,myWell.sample,myWell.detector,myWell.reporter,myWell.task,myWell.threshold,myWell.flags = tmp['well'],tmp['sample'],tmp['detector'],tmp['reporter'],tmp['task'],tmp['threshold'],dict(zip(dictKeys[17:],vals[17:]))
            res[myWell.wellNum] = myWell
        except ValueError:
            pass
        line=handle.readline()
        if not line: break 
    return res
            
    assert False, "Should not reach this line..."
        
def parseRawCycle(fname,wellData):
    """This replaces parseCycleData"""
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
        vals = np.array(map(float,vals[1:]))
        wellData[well].cycles,wellData[well].fluorData = headerVals,vals
    return

def getDetAndSamp(wellData):
    """Returns two lists of unique detectors and unique samples"""
    detectors = util.uniqify(detectors = [x.detector for x in wellData])
    samples = util.uniqify(samples = [x.sample for x in wellData])
    return detectors,samples
    
def wellIndex(data):
    index = []
    for i in range(len(data)):
        index.append(data[i]['well'])
    return index

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
    return np.log10(myArray)

#########
# Four-parameter Logistic Model fitting
#########
def nthRoot(num,n):
    return num ** (1.0/n)

def qpcrFit(self,x,a,b,x0,y0):
    """Same as fit but designed to run with optimize.curve_fit"""
    return (y0+(a/(1+((x/x0)**b))))

def qpcrFitResiduals(x,y,a,b,x0,y0):
    """
    Residuals:
    errfunc = lambda p,x,y: y-fitfunc(p,x) #Distance to the target function (residuals)
    """
    return y-qpcrFit(x,a,b,x0,y0)

def CP_FDM(p):
    return (p[2]*nthRoot(((p[1]-1)/(p[1]+1)),p[1]))

def CP_SDM(p):
    return p[2]*nthRoot((np.sqrt((3*p[1]**2)*(p[1]**2-1))-(2*(1-p[1]**2)))/((p[1]**2)+(3*p[1])+2),p[1])

def CP_SPE(p,rNoise):
    return (p[2]*nthRoot(((p[0]-rNoise)/rNoise),p[1]))

###############################
#Iterative Nonlinear Regression
###############################
def nlmFit(x,a,b,y0):
    """
    Non-linear regression function to optimize for windows in exponential phase
    here p = [a,b,y0]
    """
    return y0+(a*(b**x))

def nlmFitResiduals(x,y,a,b,y0):
    """
    Residuals:
    errfunc = lambda p,x,y: y-nlmFit(p,x) #Distance to the target function (residuals)
    """
    return y-nlmFit(x,a,b,y0)


###############################
#ddCt math
###############################
def ddCt(data,medianCts,endoControl,reference):
    tmp = {}
    #Calculate dCts
    for i in range(len(data)):
        print medianCts[data[i]['sample']]
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
    pass

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
    print vals
    vals = [i for i in vals if i != "N/A"]
    print vals
    lenvals = len(vals)
    vals.sort()
    if lenvals == 0:
        return "N/A"
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

def doPlotting(plotScript = "qPCRPlotting.q"):
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
    data = parseRawABI(mainFile)
    medianCts = aggregateReplicateCts(data) #Returns a dictionary of dictionaries by sample and then detector
    myIdx = wellIndex(data)
    
    #Efficiency Calculation from cycleFile
    print "Parsing CycleData File..."
    cycleData = parseRawCycle(cycleFile)
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
