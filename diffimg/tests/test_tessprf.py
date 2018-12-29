# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 23:11:29 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatch
#import matplotlib as mpl
#import pandas as pd
import numpy as np
import pytest

import dave.diffimg.tessprf as prf

datapath = "/home/fergal/data/tess/prf/"

def test_bracketing():
    """The method getRegularlySampledBracketingPrfs() has some fairly complicated
    bookkeeping to find the locations of the 4 prfs that brack the input col,row
    in 2d space. It internally checks this bookkeeping is correct and raises
    an assert on failure. This test exercises all 4 paths in the code.
    """
    obj = prf.TessPrf(datapath)
    
    #This raises an assertion error
    obj.getPrfAtColRow(1587, 1710, 1, 1, 1)
    obj.getPrfAtColRow(1581, 1537, 1, 1, 1)  #A
    obj.getPrfAtColRow(1579, 1537, 1, 1, 1)  #S
    obj.getPrfAtColRow(1579, 1535, 1, 1, 1)  #T
    obj.getPrfAtColRow(1581, 1535, 1, 1, 1)  #Cs


#Test out of bounds behavour
    
def test_outOfBounds():

    obj = prf.TessPrf(datapath)
    with pytest.raises(ValueError):
        obj.getPrfAtColRow(0,0, 1, 1, 1)
        
    with pytest.raises(ValueError):
        obj.getPrfAtColRow(44,1, 1, 1, 1)
        
    with pytest.raises(ValueError):
        obj.getPrfAtColRow(2093,1, 1, 1, 1)
    
    with pytest.raises(ValueError):
        obj.getPrfAtColRow(47, 0, 1, 1, 1)

    with pytest.raises(ValueError):
        obj.getPrfAtColRow(47,2048  , 1, 1, 1)
            
    #Check some in bounds
    obj.getPrfAtColRow(45, 1, 1, 1, 1)
    obj.getPrfAtColRow(2091, 2047, 1, 1, 1)
        

def testIntFloatBug():        
    """getPrfAtColRow() should return same value whether input is int or float"""

    obj = prf.TessPrf(datapath)
    
    img1 = obj.getPrfAtColRow(123, 456, 1,1,1)
    img2 = obj.getPrfAtColRow(123.0, 456.0, 1,1,1)
    
    assert np.all(img1 - img2 == 0)
    
    
def imgByOffset():
    ccd, camera, sector = 1,1,1
    col, row = 123, 456
    
    obj = prf.TessPrf(datapath)
    prfObjArray = obj.readPrfFile(ccd, camera, sector)
 
    singlePrfObj = prfObjArray[0]
    img0 = obj.getRegularlySampledPrfByOffset(singlePrfObj, 0, 0)

    for offset in range(9):
#        img1 = obj.getRegularlySampledPrfByOffset(singlePrfObj, offset, 0)
        img1 = obj.getRegularlySampledPrfByOffset(singlePrfObj, offset, 0)
        delta = img1 - img0
        
        kwargs = {'origin':'bottom', 'interpolation':'nearest', 'cmap':plt.cm.YlGnBu_r}
        plt.clf()
        plt.subplot(121)
        plt.imshow(img1, **kwargs)
        plt.colorbar()
        
        plt.subplot(122)
        kwargs['cmap'] = plt.cm.PiYG
        plt.imshow(delta, **kwargs)  
        vm = max( np.fabs( [np.min(delta), np.max(delta)] ))
#        vm = 1e-2
        plt.clim(-vm, vm)
        plt.colorbar()
        plt.suptitle(offset)
        plt.pause(.1)
        raw_input()
            
    
def testColFrac():        
    """Test that changing column fraction moves flux around"""

    obj = prf.TessPrf(datapath)
    
    img1 = obj.getPrfAtColRow(123.0, 456, 1,1,1)
    
    for frac in np.linspace(0, .9, 11):
        print("Frac is %g" %(frac))
        img2 = obj.getPrfAtColRow(123.0 + frac, 456.0, 1,1,1)
        delta = img2 - img1
        
#        prfPlot(img1, delta)
        
        #For TESS, PRFs are 13x13. Check the flux near the centre
        #is moving from lower columns to higher ones
        assert delta[6,6] >= 0, delta[6,6]
        assert delta[6,7] >= 0, delta[6,7]
        assert delta[6,5] <= 0, delta[6,5]
        

def testRowFrac():        
    """Test that changing column fraction moves flux around"""

    obj = prf.TessPrf(datapath)
    
    img1 = obj.getPrfAtColRow(123.0, 456, 1,1,1)
    
    for frac in np.linspace(0, .9, 11):
        img2 = obj.getPrfAtColRow(123.0, 456.0 + frac, 1,1,1)
        delta = img2 - img1
        
#        prfPlot(img1, delta)
        
        #For TESS, PRFs are 13x13. Check the flux near the centre
        #is moving from lower columns to higher ones
        assert delta[6,6] >= 0, delta[6,6]
        assert delta[7,6] >= 0, delta[7,6]
        assert delta[5,6] <= 0, delta[5,6]

        
def prfPlot(refImg, delta):
        
        kwargs = {'origin':'bottom', 'interpolation':'nearest', 'cmap':plt.cm.YlGnBu_r}
        plt.clf()
        plt.subplot(121)
        plt.imshow(refImg, **kwargs)
        plt.colorbar()
        
        plt.subplot(122)
        kwargs['cmap'] = plt.cm.PiYG
        plt.imshow(delta, **kwargs)  
        vm = max( np.fabs( [np.min(delta), np.max(delta)] ))
#        vm = 1e-2
        plt.clim(-vm, vm)
        plt.colorbar()
        plt.pause(.1)
    