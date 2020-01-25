# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 20:02:14 2020

Implementations of analytic PRF functions.

These are PRFs based on mathematical models of the shape of the
PRF. Two implementations available to start, a Gaussian psf,
and a Gaussian on a constant sky bacground.

"""

from ipdb import set_trace as idebug
from pdb import set_trace as debug
from scipy.special import erf
import numpy as np


#analyticprf.py
class GaussianPrf():
    """
    PRF based on a mathematical model of a Gaussian function.

    The implementation uses look up tables of Gaussian cumulative
    distribution functions for maximum speed. This is the fastest
    of all PRF functions, and a good choice when you don't
    know anything about the true shape of the stellar profile.

    This class has no configuration options.

    See also GaussianPlusSkyPrf()
    """
    #kwargs included for consistency with AbstractBaseClass
    def __init__(self, **kwargs):
        assert len(kwargs) == 0

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
    """Same as GaussianPrf, but with an additional
    constant sky background. This model seems to perform
    better for Kepler, K2 and TESS"""

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
