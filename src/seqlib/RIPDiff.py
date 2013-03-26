'''
Created on May 13, 2010

Normalizes and compares RIP vs Control (IgG or total RNA) to identify segments of transcripts that are 
preferrentially enriched in RIP

@author: lgoff
'''
##################
#Imports
##################
import intervallib
import seqstats


##################
#Classes
##################

class RIPUnit(intervallib.Interval):
    """
    Can be individual transcript or some basic unit being interrogated for differential peaks (ie. chromosome) 
    Extends intervallib.Interval class
    """
    def __init__(self,interval):
        """Initiate from existing instance of Interval class only"""
        assert isinstance(interval,intervallib.Interval)
        intervallib.Interval.__init__(interval)
        
    def scan(self):
        pass
    
    def makebins(self,binSize):
        pass
    
    def binBinom(self):
        pass
    
    def binPois(self):
        pass
    
    def fetchReads(self,bamHandle):
        pass


#################
#Functions
#################
def globalNorm(ripUnit,totReads):
    pass
    
def localNorm(ripUnitA,ripUnitB):
    pass
