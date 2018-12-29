
"""
Created on Sun Dec  2 21:03:02 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import numpy as np

import dave.diffimg.fastpsffit as ddf
import dave.diffimg.disp as disp

from dave.misc.plateau import plateau
import dave.fileio.kplrfits as kplrfits

def measurePerTransitCentroids(time, cube, period_days, epoch_days, duration_days, plotFilePattern=None):
    """Measure difference image centroids for each transit in a data set
    
    For each transit, compute the centroid of a image reprsenting the the position of the star before transit,
    after transit, and the difference image. 

    Not working yet. Before and After centroids are in the wrong place
    
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
    
#    if plotFilePattern is not None:
#        raise NotImplementedError("Plots not implemented yet")
        
        
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
            
    results = np.zeros( (len(transits), 6) )
    for i in range(len(transits)):
        key = 'transit-%04i' %(i)
        results[i,:2] = out[key]['beforeCentroid']
        results[i,2:4] = out[key]['diffCentroid']
        results[i,4:6] = out[key]['afterCentroid']
        
    out['results'] = results
    return out
    

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

    beforeGuess = pickInitialGuess(before)
    beforeSoln = ddf.fastGaussianPrfFit(before, beforeGuess)
    
    diffGuess = pickInitialGuess(diff)
    diffSoln = ddf.fastGaussianPrfFit(diff, diffGuess)

    afterGuess = pickInitialGuess(after)
    afterSoln = ddf.fastGaussianPrfFit(after, afterGuess)

    if plot:
        debuggingPlot(before, beforeSoln, beforeGuess,
                      diff, diffSoln, diffGuess,
                      after, afterSoln, afterGuess)
        #        generateDiffImgPlot(before, diff, after, beforeSoln, diffSoln, afterSoln)
#        
        plt.pause(.01)

    out = dict()
    error = 0
    error += 1 * (not beforeSoln.success)
    error += 2 * (not diffSoln.success)
    error += 4 * (not afterSoln.success)
    out['errorCode'] = error

    out['beforeCentroid'] = beforeSoln.x[:2]
    out['diffCentroid'] = diffSoln.x[:2]
    out['afterCentroid'] = afterSoln.x[:2]
    out['beforeScore'] = computeFitScore(before, beforeSoln)
    out['afterScore'] = computeFitScore(after, afterSoln)
    out['diffScore'] = computeFitScore(diff, diffSoln)


    return out


def computeFitScore(image, soln):
    
    nr, nc = image.shape
    model = ddf.computeModel(nc, nr, soln.x)
    
    score = np.sum(np.fabs( image-model ))
    score /= np.sum(image) * float(nr*nc)
    return score
        
def debuggingPlot(before, beforeSoln, beforeGuess,
                  diff, diffSoln, diffGuess,
                  after, afterSoln, afterGuess):    

    plt.clf()
    _debugPlot(before, beforeSoln, beforeGuess, 0)
    plt.subplot(3,4, 1)
    plt.ylabel("Before")
    
    _debugPlot(after, afterSoln, afterGuess, 1)
    plt.subplot(3,4, 5)
    plt.ylabel("After")

    _debugPlot(diff, diffSoln, diffGuess, 2, cmap=plt.cm.RdBu_r)
    plt.subplot(3,4, 9)
    plt.ylabel("Difference")
    
    
def _debugPlot(image, soln, guess, row, **kwargs):
    
    plt.subplot(3, 4, 4*row + 1)
    disp.plotImage(image, **kwargs)
    clim = np.min(image), np.max(image)
    plt.title("Image")
    
    plt.subplot(3, 4, 4*row + 2)
    initImage = ddf.computeModel(6,6, guess)
    disp.plotImage(initImage, clim=clim, **kwargs)
    plt.title("Initial Guess Model")
    
    plt.subplot(3, 4, 4*row + 3)
    model = ddf.computeModel(6,6, soln.x)
    disp.plotImage(model, clim=clim, **kwargs)
    plt.title("Best Fit Model")

    plt.subplot(3, 4, 4*row + 4)
    disp.plotDifferenceImage(image-model)
    plt.title("Residual")
#    plt.pause(.1)
#    raw_input()
        


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

    guess = [c0+.5, r0+.5, .5, 4*np.max(img), np.median(img)]
#    guess = [c0+.5, r0+.5, .5, 1*np.max(img), np.median(img)]
    return guess


def generateDiffImgPlot(before, diff, after, beforeSoln, diffSoln, afterSoln):
    """Generate a difference image plot"""
        
    plt.clf()
    plt.subplot(221)
    plt.title("Before")
    disp.plotImage(before)

    plt.subplot(222)
    plt.title("After")
    disp.plotImage(after)

    plt.subplot(223)
    plt.title("After - Before")
    disp.plotDiffImage(after-before)
    

    plt.subplot(224)
    plt.title("Diff")
    disp.plotDiffImage(diff)
    plotCentroidLocation(beforeSoln, 's', ms=8, label="Before")
    plotCentroidLocation(afterSoln, '^', ms=8, label="After")
    plotCentroidLocation(diffSoln, 'o', ms=12, label="Diff")
    
#    plt.legend()
    plt.pause(1)


    
        
def plotCentroidLocation(soln, *args, **kwargs):
    """Add a point to the a plot.
    
    Private function of `generateDiffImgPlot()`
    """
    col, row = soln.x[:2]
    ms = kwargs.pop('ms', 8)
    

    kwargs['color'] = 'k'
    plt.plot([col], [row], *args, ms=ms+1, **kwargs)

    color='orange'
    if soln.success:
        color='w'
    kwargs['color'] = color
    plt.plot([col], [row], *args, **kwargs)



def testSmoke():
    import dave.fileio.pyfits as pyfits
    import dave.fileio.tpf as tpf
    tic = 307210830
    sector = 2

    period_days = 3.69061
    epoch_btjd = 1356.2038
    duration_days = 1.2676/24.

    path = '/home/fergal/data/tess/hlsp_tess-data-alerts_tess_phot_%011i-s%02i_tess_v1_tp.fits'
    path = path %(tic, sector)
    fits, hdr = pyfits.getdata(path, header=True)
    
    time = fits['TIME']
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    cube = cube[:, 3:9, 2:8]
    
    res = measurePerTransitCentroids(time, cube, period_days, epoch_btjd, duration_days, "tmp")
    return res
