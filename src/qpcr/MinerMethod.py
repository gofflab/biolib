#!/usr/bin/env python

'''
Implementation of the Miner Method for qPCR crossing-point determination.

Provides the four-parameter logistic (4PL) model (``qpcrFit``), an
exponential-phase nonlinear regression model (``nlmFit``), and three
crossing-point estimation methods derived from the fitted 4PL parameters:

- FDM (First Derivative Maximum)
- SDM (Second Derivative Maximum)
- SPE (Signal-to-noise / Percentage of Efficiency)

Also contains an example fit executed at import time using a hard-coded
sample fluorescence curve (``myData``).

Reference: Zhao & Fernald (2005). "Comprehensive algorithm for
quantitative real-time polymerase chain reaction." J Comput Biol.

@author: lgoff
'''
#!/usr/bin/env python
import numpy as np

#from scipy import *
from scipy import optimize  # To do model fitting and non linear regression

# NOTE: skidmarks is not Python 3 compatible. Runs test is disabled.
# from skidmarks import wald_wolfowitz  # Required for runs test of residuals from iterative non-linear regression
#import scipy.stats.sem as sem


#myData = np.array([0.25733316,0.25389174,0.25416338,0.2587209,0.25729367,0.26071942,0.2576906,0.25828227,0.26198432,0.25957265,0.2577642,0.25586262,0.26059827,0.26065505,0.25757584,0.25949657,0.25952592,0.26461914,0.26600435,0.27098677,0.27315396,0.2857388,0.31070504,0.36050597,0.4551804,0.6308413,0.94302386,1.4290692,2.0682411,2.7252922,3.2184746,3.5508757,3.7593882,3.913022,4.034261,4.1229677,4.1557994,4.212172,4.243716,4.2849827,4.2739472,4.311232,4.322311,4.318703,4.344398])
myData = np.array([0.26943192,0.27736726,0.28434828,0.27858773,0.2779131,0.28177735,0.28615,0.2953472,0.29792145,0.30138493,0.30184093,0.30364826,0.3019202,0.3151101,0.32912096,0.34938487,0.39618066,0.4623603,0.5972733,0.84688836,1.268771,1.9334784,2.797376,3.602377,4.241921,4.687924,4.964248,5.2410073,5.3598685,5.5112166,5.6203637,5.696951,5.7454934,5.7954955,5.8482194,5.8416085,5.7862396,5.8655,5.86371,5.859713,5.874891,5.8553905,5.8210464,5.853178,5.870367])
cycles = list(map(float,range(1,len(myData)+1))) # Some platforms are fractional so I should get this from the clipped Data file.

#########
#Misc
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

#############
#Functions
############
#fitfunc = lambda p,x: p[3]+(p[0]/(1+((x/p[2])**p[1]))) # From actual paper (Zhao et al) where p = [a,b,x_0,y_0]
#errfunc = lambda p,x,y: y-fitfunc(p,x) #Distance to the target function (residuals)

def fit(p,x):
    """Evaluate the four-parameter logistic (4PL) model using a parameter vector.

    Deprecated in favor of ``qpcrFit``, which is compatible with
    ``scipy.optimize.curve_fit``.

    The model is:
        f(x) = p[3] + p[0] / (1 + (x / p[2])^p[1])

    where ``p = [a, b, x0, y0]`` following the notation in Zhao et al.

    Args:
        p: Sequence of four model parameters ``[a, b, x0, y0]``:
            a  – amplitude (difference between upper and lower asymptotes),
            b  – slope/steepness,
            x0 – inflection point (cycle at midpoint),
            y0 – baseline fluorescence (lower asymptote).
        x: Cycle number (scalar or array).

    Returns:
        Predicted fluorescence value(s) at cycle ``x``.
    """
    return (p[3]+(p[0]/(1+((x/p[2])**p[1]))))

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
for i in range(len(p0)):
    pSEC.append(np.sqrt(pCov[i][i]))

#RNoise is standard error of y0
RNoise = pSEC[3]

print(p0)
print(p1)
print(RNoise)
print(CP_FDM(p1))
print(CP_SDM(p1))
print(CP_SPE(p1,RNoise))
#print(myData)
#print(fitData)
print("###############")

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
# NOTE: runsTest is disabled because skidmarks is not Python 3 compatible.
# runsTest = wald_wolfowitz(run)
pass  # runsTest disabled

print(lmParams)
print(xdata)
print(ydata)
print(lmFitData)
print(lmResids)

print("#################")
print(run)
# print(1-runsTest['p'])  # runsTest disabled
