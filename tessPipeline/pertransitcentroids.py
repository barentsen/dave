
"""
Created on Sun Dec  2 21:03:02 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl
import pandas as pd
import numpy as np

import dave.diffimg.fastpsffit as ddf
from kepler.plateau import plateau
import kplrfits


def measurePerTransitCentroids(time, cube, period_days, epoch_days, duration_days, plotFilePattern=None):
    """Measure difference image centroids for each transit in a data set
    
    For each transit, compute the centroid of a image reprsenting the the position of the star before transit,
    after transit, and the difference image. 
    
    Inputs
    --------
    time
        (1d np array) Times of cadences
    cube
        (3d np array) Contents of TPF file for this target
    period_days
        (float) Period of transit
    epoch_days
        (float) Time of first transit
    duration_days
        (float) Duration of transit. Results aren't strongly dependent on whether you choose total
        transit time, or just from end of ingress to start of egress, but the latter is preferable
        
    Optional Inputs
    -------------
    plotFilePattern
        If **None**, genereate no diagnostic plots. Otherwise, save plots to this path and file root
        
    Returns
    ---------
    A dictionary.
    """
    plot = True
    if plotFilePattern is None:
        plot = False
    
    if plotFilePattern is not None:
        raise NotImplementedError("Plots not implemented yet")
        
        
    isnan = np.isnan(time)
    time = time[~isnan]
    cube = cube[~isnan]

    out = dict()
    transits = getIngressEgressCadences(time, period_days, epoch_days, duration_days)
    for i in range(len(transits)):
        cin = transits[i]
        res = measureCentroidShift(cube, cin, plot)
        out['transit-%04i' %(i)] = res
        

        if plotFilePattern is not None:
            plt.suptitle('%s-trans%02i' %(plotFilePattern , i))   
            plt.savefig('%s-trans%02i.png' %(plotFilePattern , i))

    results = np.array( (len(transits), 6) )
    for i in range(len(transits)):
        key = 'transit-%04i' %(i)
        results[:,2] = out[key]['beforeCentroid']
        results[:,2:4] = out[key]['diffCentroid']
        results[:,4:6] = out[key]['afterCentroid']
        
    out['results'] = results
    
    

def getIngressEgressCadences(time, period_days, epoch_btjd, duration_days):
    """Get a list of transit start and end times in units of cadence number
    
    Inputs
    ----------
    
    
    Returns
    ----------
    A 2d numpy array. zeroth column is cadence number of transit starts, first column is cadence
    number of transit ends.
    """
    assert np.all(np.isfinite(time))

    idx = kplrfits.markTransitCadences(time, period_days, epoch_btjd, duration_days)
    transits = np.array(plateau(idx, .5))

    return transits


def measureCentroidShift(cube, cin, plot=True):
    """
    
    
    Todo
    ---------
    * Adapt this code so it can be parameterised to accept different fitting algorithms, instead
      of hard coding fastGaussianPrf
    """
    
    before, after, diff = generateDiffImg(cube, cin, plot=plot)

    print("Before...")
    guess = pickInitialGuess(before)
    beforeSoln = ddf.fastGaussianPrfFit(before, guess)

    print("Diff...")
    guess = pickInitialGuess(diff)
    diffSoln = ddf.fastGaussianPrfFit(diff, guess)

    print("After...")
    guess = pickInitialGuess(after)
    afterSoln = ddf.fastGaussianPrfFit(after, guess)


    out = dict()
    error = 0
    error += 1 * (not beforeSoln.success)
    error += 2 * (not diffSoln.success)
    error += 4 * (not afterSoln.success)
    out['errorCode'] = error

    out['beforeCentroid'] = beforeSoln.x[:2]
    out['diffCentroid'] = diffSoln.x[:2]
    out['afterCentroid'] = afterSoln.x[:2]

    return out


def generateDiffImg(cube, transits, offset_cadences=3, plot=False):
    """Generate a difference image.

    Also generates an image for each the $n$ cadedences before and after the transit,
    where $n$ is the number of cadences of the transit itself

    Inputs
    ------------
    cube
        (np 3 array) Datacube of postage stamps
    transits
        (2-tuples) Indices of the first and last cadence

    Optional Inputs
    -----------------
    offset_cadences
        (int) How many cadences gap should be inserted between transit and before and after images 
        Bryson et al. suggest 3 is a good number.
    plot
        (Bool) If true, generate a diagnostic plot


    Returns
    -------------
    Three 2d images,

    before
        The sum of the n cadences before transit (where n is the number of in-transit cadences
    after
        The sum of the n cadences after transit
    diff
        The difference between the flux in-transit and the average of the flux before and after

    Notes
    ---------
    * When there is image motion, the before and after images won't be identical, and the difference
    image will show distinct departures from the ideal prf.
    
    Todo
    ---------
    * This function belongs more properly in diffimg.py. I should port it there.
    """

    dur  = transits[1] - transits[0]
    s0, s1 = transits - dur - offset_cadences
    e0, e1 = transits + dur + offset_cadences

    before = cube[s0:s1].sum(axis=0)
    during = cube[transits[0]:transits[1]].sum(axis=0)
    after = cube[e0:e1].sum(axis=0)

    diff = .5 * (before + after) - during
    return before, after, diff    


def pickInitialGuess(img):
    """Pick initial guess of params for `fastGaussianPrfFit`
    
    Inputs
    ---------
    img
        (2d np array) Image to be fit
        
    Returns
    ---------
    An array of initial conditions for the fit
    """
    r0, c0 = np.unravel_index( np.argmax(img), img.shape)

    guess = [c0+.5, r0+.5, .5, np.max(img), np.median(img)]
    return guess



def generateDiffImgPlot(before, diff, after):
    """Generate a difference image plot"""
    plt.clf()
    plt.subplot(221)
    plt.imshow(before, origin='bottom')
    plt.title("Before")
    plt.colorbar()

    plt.subplot(222)
    plt.imshow(after, origin='bottom')
    plt.title("After")
    plt.colorbar()

    plt.subplot(223)
    plt.imshow(after - before, origin='bottom', cmap=plt.cm.RdYlBu_r)
    plt.title("After - Before")
    plt.colorbar()

    plt.subplot(224)
    plt.imshow(diff, origin='bottom', cmap=plt.cm.RdYlBu_r)
    plt.title("Diff")
    plt.colorbar()


#def main(tic, sector, period_days, epoch_btjd, duration_days, outpattern):
#
#    path = '/home/fergal/data/tess/hlsp_tess-data-alerts_tess_phot_%011i-s%02i_tess_v1_tp.fits'
#    path = path %(tic, sector)
#    fits, hdr = pyfits.getdata(path, header=True)
#    cube = ktpf.getTargetPixelArrayFromFits(fits, hdr)
#    cube = cube[:, 3:9, 2:8]
#
#    time = fits['TIME']
#    isnan = np.isnan(time)
#    time = time[~isnan]
#    cube = cube[~isnan]
#
#    transits = getIngressEgressCadences(time, period_days, epoch_btjd, duration_days)
#    with open('%s.cent.txt' %(outpattern), 'w') as fp:
#        for i in range(len(transits)):
#            print("Transit %i" %(i))
#            cin = transits[i]
#            res = measureCentroidShift(cube, cin, True)
#
#            plt.suptitle('%s-trans%02i' %(outpattern, i))
#            plt.savefig('%s-trans%02i.png' %(outpattern, i))
#
#            pattern = "%.6f " * len(res)
#            pattern = pattern + "\n"
#            fp.write( pattern % tuple(res))
