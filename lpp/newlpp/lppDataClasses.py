#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:14:04 2018

@author: smullally
"""
from __future__ import print_function
from __future__ import division
import scipy.io as spio
from astropy.io import fits
import requests
import numpy as np


class TCE(object):
    
    def __init__(self, catid, planetNum):
        """
        Init creates default values that define the TCE.
        The minimum set are (id, times (days), fluxes (fractional), period (days), tzero)
        """
        self.id = catid
        self.planetNum = planetNum
        self.time = np.array([0])
        self.phase = np.array([0])
        self.flux = np.array([0])
        
        self.period = -1  #period in days
        self.tzero = 0
        self.dur = -1  #duration in hours
        self.depth = 0  #transit depth in ppm.
        self.mes = 10
        
    def populateFromDvExt(self, data, header):
        """
        Fill in the TCE informaiton using a DVT data extension and header.
        This uses the median detrended light curve.t
        """
        
        self.time = data['time']
        self.phase = data['phase']
        self.flux = data['lc_detrend']
        self.period = header['TPERIOD']
        self.tzero = header['TEPOCH']
        self.dur = header['TDUR']
        self.depth = header['TDEPTH']
        self.mes = header['maxmes']
        
        self.checkTce()
        
    def checkTce(self):
        """
        Check basic properties of the TCE values to ensure valid.
        """
        
        if len(self.time) != len(self.phase):
            raise Warning("Length of Time and Phase do not agree.")
            
        if len(self.time) != len(self.flux):
            raise Warning("Length of Times and Fluxes do not agree.")
        
        if self.period <= 0 :
            raise Warning("Period in days has a value of zero or less.")
        
        if self.dur <= 0 :
            raise Warning("Duration in hrs has a value of zero or less.")


#----

class MapInfo(object):
    """
    Class to read in the matlab blob with the map.
    """
    
    def __init__(self,filename):
        
        self.filename=filename
        
        self.readMatlabBlob(filename)

    
    def readMatlabBlob(self,filename):
        """
        read in matlab blob
        Using the DV trained one.
        """      

        mat=spio.loadmat(filename,matlab_compatible=True)
        
        #Pull out the information we need.
        
        self.n_dim = mat['mapInfoDV']['nDim'][0][0][0][0]
        self.Ymap = mat['mapInfoDV']['Ymap'][0][0][0][0]
        self.YmapMapping = self.Ymap['mapping']
        self.YmapMean = self.YmapMapping['mean'][0][0][0]
        self.YmapM = self.YmapMapping['M'][0][0]
        self.YmapMapped = self.Ymap['mapped']
        self.knn=mat['mapInfoDV']['knn'][0][0][0][0]
        self.knnGood=mat['mapInfoDV']['knnGood'][0][0][:,0]
        self.mappedPeriods=mat['mapInfoDV']['periods'][0][0][0]
        self.mappedMes=mat['mapInfoDV']['mes'][0][0][0]
        self.nPsample=mat['mapInfoDV']['nPsample'][0][0][0][0]  #number to sample
        self.nPercentil=mat['mapInfoDV']['npercentilTM'][0][0][0][0]
        self.dymeans=mat['mapInfoDV']['dymean'][0][0][0]
        self.ntrfr= 2.0
        self.npts=80.0
        
        
