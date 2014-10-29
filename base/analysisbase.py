#!/usr/bin/env python
"""
*******************************************************
AnalysisBase for histogram-based analysis.

pascal nef, Oct 28th 2014
*******************************************************
"""
import sys
import os
import re
import math
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


class AnalysisBase(object):
    """ Base class for analysis """
    def __init__(self):
        self.H = {}
        self.fig = 0
        self.figID=0

    def gethistofromarray(self, ID, data, **kwargs): 
        """ create np.histogram from array """

        try:
            self.H[ID]=Histo( ID)
            self.H[ID].gethistofromarray(data, **kwargs)
        except:
            print "could not create histo with kwargs:"
            for name, value in kwargs.items():
                print name, value
            raise

    def makelikelihood(self, ID1, ID2, IDnew):
        """ create likelihood based on signal (ID1) and background (ID2) histos"""
        # check if binning of two histos is identical
        try:
            (self.H[ID1].bins == self.H[ID2].bins).all()
        except:
            print "histos have different binning:"
            print ID1, self.H[ID1].bins 
            print ID2, self.H[ID2].bins 
        
        self.H[IDnew]=Histo(IDnew)
        self.H[IDnew].hist = np.zeros(len(self.H[ID1].hist ), dtype=np.float)
        self.H[IDnew].bins = self.H[ID1].bins
        for i,(s,b)  in enumerate(zip(self.H[ID1].hist, (self.H[ID2].hist))):
            if (s+b)>0:
                self.H[IDnew].hist[i]=float(s)/float(s+b)



    def evaluate_likelihood(self, data, likelihoods):
        """ evaluate the data for a given (list of) likelihoods
            -> dimension of the data must correspond to the number of likelihoods
            -> i.e. if data = [(a1,b1), (a2,b2)], result is [L1(a1)*L2(b1), L1(a2)*L2(b2)
            -> this is a projective likelihood in case dimension >1
        """    
        evaluated=[]
        for point in data:
            result = 1
            out_of_range = False
            for dim, value in enumerate(point):
                likelihood_value = self.evaluate_hist_at_point(value, likelihoods[dim])
                if(likelihood_value==-1):
                    out_of_range = True
                else:
                    result *= likelihood_value
            if not out_of_range:
                evaluated.append(result)
                    
        return evaluated
                


    def evaluate_hist_at_point(self, point, ID):
        """ evaluate a histogram at a certain point along x-axis"""
        if(point < self.H[ID].bins[0] or point >= self.H[ID].bins[-1]):
            return -1
        else:
            binnumber = self.H[ID].bins.searchsorted(point, side='right')-1
            return self.H[ID].hist[binnumber]
    



class Histo(object):
    """ Histogram base class """
    def __init__(self, histname):
        self.histname       = histname
        self.markercolor    = "b"
        self.linecolor      = "b"
        self.markerstyle    = "o"
        self.hist           = []
        self.bins           = []

    def gethistofromarray(self, data, **kwargs):
        # parse args -----------
        for key in ('range', 'bins', 'density'):
            if key in kwargs:
                setattr(self, key, kwargs[key])

        self.hist ,self.bins = np.histogram(data, **kwargs)




def wait():
    var = raw_input("waiting. hit enter to continue")
