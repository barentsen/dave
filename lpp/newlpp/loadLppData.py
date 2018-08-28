#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 11:58:23 2018

Classes to read in data files containing either the 
light curve containing the transit (Timeseries) or
the information about the mapping (MapInfo)
These are used by lpp_transform

@author: smullally
"""

import scipy.io as spio
from astropy.io import fits


class TCE(object):
    
    def __init__(self, filename, ext=1):
        self.filename = filename
        self.ext=ext
                
    def readDV(self):

        try:
            hdu=fits.open(self.filename)
        except IOError:
            print "Filename not found."
            
        ext=self.ext
        
        self.time=hdu[ext].data['TIME']
        self.flux=hdu[ext].data['LC_DETREND']
        self.period=hdu[ext].header['TPERIOD']
        self.tzero=hdu[ext].header['TEPOCH']
        self.dur=hdu[ext].header['TDUR']
        self.depth=hdu[ext].header['TDEPTH']
        self.mes=hdu[ext].header['MAXMES']
        
        hdu.close()

class MapInfo(object):
    
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
        self.YmapMapping = self.Ymap['mapping'][0][0][0]
        self.YmapMapped = self.Ymap['mapped']
        self.knn=mat['mapInfoDV']['knn'][0][0][0][0]
        self.knnGood=mat['mapInfoDV']['knnGood'][0][0][:,0]
        self.mappedPeriods=mat['mapInfoDV']['periods'][0][0][0]
        self.nPsample=mat['mapInfoDV']['nPsample'][0][0][0][0]  #number to sample
        self.nPercentil=mat['mapInfoDV']['npercentilTM'][0][0][0][0]
        self.dymeans=mat['mapInfoDV']['dymean'][0][0][0]
        self.ntrfr= 2.0
        self.npts=80.0
