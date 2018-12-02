# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Fri Nov 30 23:11:29 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpatch
#import matplotlib as mpl
#import pandas as pd
import numpy as np
import pytest

import prf

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
        
    
def testColFrac():        
    """Test that changing column fraction moves flux around"""
    assert False