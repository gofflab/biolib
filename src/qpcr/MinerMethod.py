#!/usr/bin/env python

'''
Created on Sep 1, 2010

@author: lgoff
'''
#!/usr/bin/env python
import numpy as np
#from scipy import *
from scipy import optimize # To do model fitting and non linear regression
from skidmarks import wald_wolfowitz # Required for runs test of residuals from iterative non-linear regression
#import scipy.stats.sem as sem


#myData = np.array([0.25733316,0.25389174,0.25416338,0.2587209,0.25729367,0.26071942,0.2576906,0.25828227,0.26198432,0.25957265,0.2577642,0.25586262,0.26059827,0.26065505,0.25757584,0.25949657,0.25952592,0.26461914,0.26600435,0.27098677,0.27315396,0.2857388,0.31070504,0.36050597,0.4551804,0.6308413,0.94302386,1.4290692,2.0682411,2.7252922,3.2184746,3.5508757,3.7593882,3.913022,4.034261,4.1229677,4.1557994,4.212172,4.243716,4.2849827,4.2739472,4.311232,4.322311,4.318703,4.344398])
myData = np.array([0.26943192,0.27736726,0.28434828,0.27858773,0.2779131,0.28177735,0.28615,0.2953472,0.29792145,0.30138493,0.30184093,0.30364826,0.3019202,0.3151101,0.32912096,0.34938487,0.39618066,0.4623603,0.5972733,0.84688836,1.268771,1.9334784,2.797376,3.602377,4.241921,4.687924,4.964248,5.2410073,5.3598685,5.5112166,5.6203637,5.696951,5.7454934,5.7954955,5.8482194,5.8416085,5.7862396,5.8655,5.86371,5.859713,5.874891,5.8553905,5.8210464,5.853178,5.870367])
cycles = map(float,range(1,len(myData)+1)) # Some platforms are fractional so I should get this from the clipped Data file.

#########
#Misc
#########
def nthRoot(num,n):
    return num ** (1.0/n)

#############
#Functions
############
#fitfunc = lambda p,x: p[3]+(p[0]/(1+((x/p[2])**p[1]))) # From actual paper (Zhao et al) where p = [a,b,x_0,y_0]
#errfunc = lambda p,x,y: y-fitfunc(p,x) #Distance to the target function (residuals)

def fit(p,x):
    """
    Depricated in favor of qpcrFit to use optimize.curve_fit()
    f(x) Logistic model for qPCR Data
    fitfunc = lambda p,x: p[3]+(p[0]/(1+((x/p[2])**p[1]))) # From actual paper (Zhao et al) where p = [a,b,x_0,y_0]
    """
    return (p[3]+(p[0]/(1+((x/p[2])**p[1]))))

def qpcrFit(x,a,b,x0,y0):
    """Same as fit but designed to run with optimize.curve_fit"""
    return (y0+(a/(1+((x/x0)**b))))

def qpcrFitResiduals(x,y,a,b,x0,y0):
    """
    Residuals:
    errfunc = lambda p,x,y: y-fitfunc(p,x) #Distance to the target function (residuals)
    """
    return y-qpcrFit(x,a,b,x0,y0)

def nlmFit(x,a,b,y0):
    """
    Non-linear regression function to optimize for windows in exponential phase
    here p = [a,b,y0]
    """
    return y0+(a*(b**x))

def nlmFitResiduals(x,y,a,b,y0):
    """
    Residuals:
    errfunc = lambda p,x,y: y-nlmFit(x,a,b,y0) #Distance to the target function (residuals)
    """
    return y-nlmFit(x,a,b,y0)

def CP_FDM(p):
    return (p[2]*nthRoot(((p[1]-1)/(p[1]+1)),p[1]))

def CP_SDM(p):
    return p[2]*nthRoot((np.sqrt((3*p[1]**2)*(p[1]**2-1))-(2*(1-p[1]**2)))/((p[1]**2)+(3*p[1])+2),p[1])

def CP_SPE(p,rNoise):
    return (p[2]*nthRoot(((p[0]-rNoise)/rNoise),p[1]))


Y0 = np.mean(myData[:5]) # Initial guess as to baseline fluorescence (mean of first five cycles)
X0 = cycles[np.argmin(abs(myData-np.mean(myData)))] # Initial guess as to inflection point at middle of curve
a = (np.max(myData)-np.min(myData)) # Initial guess as to y value at inflection of 
b = 0

#p0 = [np.mean(myData[:5]),2.,median(myData),np.mean(myData[-5:])]
p0 = [a,b,X0,Y0]

#p1,success = optimize.leastsq(residuals,p0,args = (cycles,myData),maxfev=5000)
#p1,pCov = myOptimize(func=qpcrFit,xdata=cycles,ydata=myData,p=[a,b,X0,Y0])
p1,pCov = optimize.curve_fit(qpcrFit,xdata=cycles,ydata=myData,maxfev=5000)

fitData = [qpcrFit(x,p1[0],p1[1],p1[2],p1[3]) for x in cycles]

pSEC = []
#Get standard error of regression coefficients
for i in xrange(len(p0)):
    pSEC.append(np.sqrt(pCov[i][i]))

#RNoise is standard error of y0
RNoise = pSEC[3]

print p0
print p1
print RNoise
print CP_FDM(p1)
print CP_SDM(p1)
print CP_SPE(p1,RNoise)
#print myData
#print fitData
print "###############"

#Iterative Nonlinear Regression
i = 15
j = 21
xdata=range(i,j)
ydata=myData[i+1:j+1]

lmParams,lmCov = optimize.curve_fit(nlmFit,xdata=xdata,ydata=ydata,maxfev=5000)
lmFitData = [nlmFit(x,lmParams[0],lmParams[1],lmParams[2]) for x in xdata]
lmResids = nlmFitResiduals(xdata,ydata,lmParams[0],lmParams[1],lmParams[2])

#P-value for runs test on resids
run = [x>=0 for x in lmResids]
runsTest = wald_wolfowitz(run)

print lmParams
print xdata
print ydata
print lmFitData
print lmResids

print "#################"
print run
print 1-runsTest['p']