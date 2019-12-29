# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 11:06:01 2019

@author: fergal
"""

from pdb import set_trace as debug
import numpy as np


def markTransitCadences(time, period_days, epoch_bkjd, duration_days,\
    numberOfDurations=1, flags=None):
    """Create a logical array indicating which cadences are
    affected by a transit

    Input:
    -------------
    time:
        (numpy 1d array) array of cadence times
    period_days:
        Transit period
    epoch_bkjd:
        Transit epoch
    duration_days:
        Duration of transit (start to finish). If you
        select a duration from first to last contact,
        all cadences affected by transit are selected.
        If you only select 2 to 3 contact, only the
        interior region of the transit is selected
    numberOfDurations
        How much of the lightcurve either side of
        the transit to mark. Default means to mark
        1/2 the transit duration either side of
        transit center.

    Optional Inputs:
    flags:
        If set, must be an array of booleans of length ``time``.
        Cadences where this array is true are ignored in the
        calculation. This is useful if some of the entries of time
        are Nan.

    Returns:
    -------------
    Array of length len(time) of booleans. Element set to true
    are affected by transits
    """

    if flags is None:
        flags = np.zeros_like(time, dtype=bool)

    i0 = np.floor((np.min(time[~flags])-epoch_bkjd)/period_days)
    i1 =  np.ceil((np.max(time[~flags])-epoch_bkjd)/period_days)
    assert(np.isfinite(i0))
    assert(np.isfinite(i1))

    irange = np.arange(i0, i1+1)
    transitTimes = epoch_bkjd + period_days*irange

    idx = np.zeros( len(time), dtype=np.bool8)
    for tt in transitTimes:
        diff = time - tt
        diff[flags] = 1e99  #A large value that isn't Nan
        assert np.all(np.isfinite(diff)), "Nan found in diff"

        idx = np.bitwise_or(idx, np.fabs(diff) < \
            .5*duration_days*numberOfDurations)

    if not np.any(idx):
        print ("WARN: No cadences found matching transit locations")
    return idx



def medianZero(data, col):
    """If x= data[:,col], set x = (x/median(x) - 1

        Note:
        This function is preferable to meanZero when the data has outliers.
    """

    med = np.median(data[:,col])

    tol=1e-9
    if np.fabs(med) < tol:
        raise ValueError("Median of lightcurve already nearly zero")

    data[:,col] /= med
    data[:,col] -= 1
    return data


def getIndicesOfBadData(time, flux, flags, mask):
    return ~getIndicesOfGoodData(time, flux, flags, mask)


def getIndicesOfGoodData(time, flux, flags, mask):
    """Remove data that is flagged as bad or has nans for time or flux

    Kepler data can be bad in three ways. The flux can be set to NaN,
    the *time* can be set to NaN, or a flag can be set. Not
    all flags are equally bad, so we include data flagged as
    zero crossing, spsd, or cosmic rays in the calibration data.
    These data may indicate lower quality data.


    Input:
    time
        (1d np array)
    flux
        (1d np array)
    flags
        (1d np array)
    mask
        (int)
        
    Returns
    An array of boolean containing the indices of good data
    """

    assert len(time) == len(flux)
    assert len(time) == len(flags)
    
    flags = flags.astype(np.uint32)
    
    #Set a mask that matches everything but listed flags.
    mask = np.zeros( (len(flags)), dtype=np.uint32) + mask

    #Find any cadences that match flag
    idx = flags & mask
    idx = flags > 0

    #Also match cadences where we don't flag but set either
    #time or flux to NaN
    idx &= np.isfinite(time)
    idx &= np.isfinite(flux)

    #Return only the indices set to true
    return idx

