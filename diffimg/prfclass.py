# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:43:43 2020

A sketch at a unified API for accessing the various PRF/PSF classes


The idea should be that every daughter class implements getPrfAtColRow(),
and the abstract class implements getPrfForBoundingBox.

The abstract lookup class will implement any shared bookkeeping code
for interpolation for PRFs.


File layout

abstractprf.py
    AbstractPrf
    AbstractLookupPrf

kepler.py
    KeplerPrf
    K2Prf  #Mneumonic?

tess.py
    TessPrf

analyticprf.py
    GaussianPrf
    GaussianPlusSkyPrf
    GeneralPrf  #For fitting arbitrary model

centroid.py
    fitPrfToImage
    fitPrfToCube
    measureDiffImageOffset

k2centroid
    #Diff images for K2 are special because of the roll
    measureK2DiffImageOffset

@author: fergal
"""

from ipdb import set_trace as idebug
from pdb import set_trace as debug
import numpy as np

from scipy.special import erf


#abstractprf.py
class AbstractPrf():
    def __init__(**kwargs):
        pass

    def getPrfForBbox(self, bbox, col0, row0, **kwargs):
        pass


class AbstractLookupPrf():
    """Abstract class for PRFs that are computed by lookup of the data in a file"""

    def _interpolatePrf(self, bbox, col, row, getPrfFunc):
        """Given 4 PRF objects at 4 locations on CCD, figure out
        what the PRF at our location should look like

        This method will be called by any daughter class from withing
        `getPrfForBbox`
        """
        raise NotImplementedError



#analyticprf.py
class GaussianPrf():
    """See also GaussianPlusSkyPrf()"""
    def __init__(self, **kwargs):
        assert len(kwargs) == 0
        pass

    def getPrfForBbox(self, bbox, col0, row0, sigma=1, flux=1):
        c0, c1, r0, r1 = bbox
        nc = c1 - c0
        nr = r1 - r0

        xc = np.arange(nc) + c0
        xr = np.arange(nr) + r0
        cols, rows = np.meshgrid(xc, xr)

        model = _analytic_gaussian_integral(cols, rows, col0, row0, sigma, flux)
        return model


class GaussianPlusSkyPrf(GaussianPrf):
    """An analytic PRF modeling the star as a 2d symettric Gaussian
    sitting on top of a constant background"""

    def getPrfForBbox(self, bbox, col0, row0, sigma=1, flux=1, sky=0):
        model = GaussianPrf.getPrfForBbox(bbox, col0, row0, sigma, flux)
        return model + sky


def _analytic_gaussian_integral(self, col, row, col0, row0, sigma0, flux0):
    """Compute the flux per pixel for a symetric Gaussian centred at
    col0, row0
    """
    z_col1 = .5 * (col - col0) / sigma0
    z_col2 = .5 * (col+1 - col0) / sigma0

    z_row1 = .5 * (row - row0) / sigma0
    z_row2 = .5 * (row+1 - row0) / sigma0

    flux = flux0
    flux *= phi(z_col2) - phi(z_col1)
    flux *= phi(z_row2) - phi(z_row1)
    return flux


sqrt2 = np.sqrt(2)
def phi(z):
    """Compute integral of gaussian function in the range (-Inf, z],

    `z` is defined as (x - x0) / sigma, where x0 is the central value of the Gaussian.

    See `scipy.special.erf` for details
    """
    return .5 * ( 1 + erf(z/sqrt2) )



#keplerprf.py
class KeplerPrf(AbstractLookupPrf):
    def __init__(self, path, **kwargs):
        AbstractLookupPrf.__init__(self, path)
        self.module = kwargs['module']
        self.output = kwargs['output']
        self.gridSize = 50

    def getPrfForBbox(self, bbox, col0, row0, **kwargs):
        pass


class K2Prf(KeplerPrf):
    """Currently identical to the KeplerPrf, but I include this class
    in case that changes in the future"""
    pass


#tessprf.py
class TessPrf(AbstractLookupPrf):
    def __init__(self, path, **kwargs):
        AbstractLookupPrf.__init__(self, path)
        self.sector = kwargs['sector']
        self.camera = kwargs['camera']
        self.ccd = kwargs['ccd']
        self.gridSize = 9

    def getPrfForBbox(self, bbox, col0, row0, sector, camera, ccd):
        pass


#centroid.py

def fitPrfToCube(prfObj, cube, bbox, *args):
    """Find the best fit prf for each image in cube.

    See fitPrfToImage for more details

    args is [col, row, [param1, [param2, ...]]]


    Notes
    ------
    getPrfForBBox(bbox, col0, row0, *args)
    costFun could argubly be an input parameter
    """

    assert len(args) >= 2

    initGuess = args
    options = {'disp':False, 'eps':.02, 'maxiter':80}
    #Don't let the fit leave the bounding box
    bounds=[(bbox[0], bbox[1]), \
            (bbox[2], bbox[3]), \
            (1, None), \
           ]

    #Define bounds for other params

    solnArray = []
    for i in range(len(cube)):
        img = cube[i]
        argsForCostFunc = (prfObj, img, bbox)

        #@TODO update initial guess if previous fit succeeded


        import scipy.optimize as spopt
        #This doesn't work, can't pass dict to minimize
        soln = spopt.minimize(costFunc, initGuess, argsForCostFunc, \
            method="L-BFGS-B", bounds=bounds, options=options)

        solnArray.append(soln)

    return solnArray


def costFunc(args, prfObj, img, bbox):
    """Measure goodness of fit (chi-square) for a PRF
    with a given col, row and scale (i.e brightness)


    args is [col, row, [param1, [param2, ...]]]
    """

    col = x[0] - bbox[0]
    row = x[1] - bbox[2]  #Note, two, not one
    model = prfObj.getPrfForBbox(bbox, *args)

    cost = img-model
    cost = np.sum(cost**2)
    return cost


def fitPrfToImage(prfObj, img, bbox, col0, row0, params):
    """Find the best fit prf for a star centred at col0, row0 in image

    Typically, the PRF model is fit to a single star in a larger image.
    It makes little sense to fit the entire image, when only a single star
    is of interest. Instead we cut out a "postage stamp", or a sub image
    and fit that. The bounding box `bbox` is a tuple of 4 integers representing
    the pixel coordinates of the edges of the postage stamp. For example,
    if the sub images is created by selecting::

        subImg = fullImage[100:120, 300:320]
        bbox = (300, 320, 100, 120)  #(c1, c2, r1, r2)

    Note how the numpy array is addressed as row column, but the bounding box
    is constructed as columns, then rows. If you are in fact fitting the
    entire image, set bbox to::

        nr, nc = img.shape
        bbox = (0, nc, 0, nr)


    Inputs
    ----------
    prfObj
        A daughter class of AbstractPrf, or one that implements a
        `getPrfForBbox()` method

    img
        (2d numpy array) Array of values representing an image to fit

    bbox
        (tuple of 4 integers) Pixel coordinates of the boundararies of
        the image. See discussion above

    col0, row0
        (floats) Initial guess as to the centroid of the star. Typically,
        these values should be inside the bounding box, but aren't required
        to be. It is sometimes unavoidable that the centroid is outside the
        bounding box

    params
        (dictionary) Other fittable parameters. Keys and values depend on
        the type of PRF object being fitted. For example, the Gaussian PRF
        takes arguments of sigma, background level, etc, while the Kepler
        PRF takes no other tunable params.

    Returns
    ---------
    A scipy.optimise.OptimizeResult object. The best fit
    params are stored in OptimizeResult.x

    """

    #TODO: Do checks on inputs

    cube = np.atleast_3d(img)
    solnArray = fitPrfToCube(prfObj, cube, col0, row0, params)
    return solnArray[0]


