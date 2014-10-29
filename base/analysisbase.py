#!/usr/bin/env python
##################################################################
# pascal nef                                July 25, 2014        #
##################################################################
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
    def __init__(self):
        self.H = {}
        self.fig = 0
        self.figID=0

    def gethistofromarray(self, ID, data, **kwargs): 

        try:
            self.H[ID]=Histo( ID)
            self.H[ID].gethistofromarray(data, **kwargs)
        except:
            print "could not create histo with kwargs:"
            for name, value in kwargs.items():
                print name, value
            raise

    def makelikelihood(self, ID1, ID2, IDnew):
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
        if(point < self.H[ID].bins[0] or point >= self.H[ID].bins[-1]):
            return -1
        else:
            binnumber = self.H[ID].bins.searchsorted(point, side='right')-1
            return self.H[ID].hist[binnumber]
    


    def NewFig(self):
        self.fig= plt.figure(num=self.figID, figsize=(4, 4), dpi=200)
        self.figID +=1
    

class Histo(object):
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


    def setp(self, **kwargs):
        try: 
            self.markerstyle = kwargs['marker']
        except:
            pass
        try:
            self.markercolor = kwargs['markercolor']
        except:
            pass
        try:
            self.linecolor = kwargs['linecolor']
        except:
            pass




def wait():
    var = raw_input("waiting. hit enter to continue")
