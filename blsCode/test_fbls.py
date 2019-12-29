# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:21:58 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import fbls2 as fbls

def test_smoke():
    """Does fbls crash when given sane inputs?"""

    size = 1000
    t = np.arange(size)/48.  #48 cadences per day
    y = np.random.randn(size)

    blsObj = fbls.BlsSearch(t, y, 1, [1, 10], [2/24., 4/24.])

    result = blsObj.getEvent()
    assert(len(result) == 5)



def test_strongTransit():

    size = 1000
    t = np.arange(size)/48.  #48 cadences per day
    y = np.random.randn(size)

    period = 5.5  #days
    t0 = 2.3
    depth = 5  #% sigma signficant depth per *point*. A v strong transit

    idx = np.fabs(np.fmod(t+t0, period)) < .25
    y[idx] -= depth

    mp.clf()
#    mp.plot(t, np.fmod(t +t0, period), 'k-')
#    mp.plot(t, y, 'k.')
#    mp.plot(t[idx], y[idx], 'ro')

    blsObj = fbls.BlsSearch(t, y, 1, [1, 10], [4/24., 6/24., 10/24.])
    result = blsObj.getEvent()

    blsObj.plot()
    print result



import mpmath as mpm

def computeLookupTable():
    mpm.mp.dps = 1000

    for z in range(100):
        print "%i %.2f" %(z, mpm.log10( mpm.log10( mpm.ncdf(z) )))