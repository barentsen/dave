
"""
Created on Fri Feb 15 15:27:38 2019

Compute the chi-squared of the BLS transit fit to the data
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

from dave.fileio.kplrfits import markTransitCadences
import dave.misc.noise as noise


def computeBlsChiSq(time_days, flux_frac, period_days, epoch_days, duration_hrs, depth_ppm):

    duration_days = duration_hrs / 24.

    return computeBlsChiSqUnitless(time_days, flux_frac,
                                   period_days, epoch_days,
                                   duration_days, depth_ppm)


def computeBlsChiSqUnitless(time, flux, period, epoch, duration, depth):
    """Compute the reduced chi squared of the box transit fit to the data

    This function is a simplified version of the chi-squared veto used by
    the Kepler SOC pipeline. The idea is that if the BLS is fooled by a
    couple of outlier points, the chi-squared of the fit of the box transit
    to the folded data will be extremely high.

    The SOC version requires a whole paper to describe the method, we think
    we can get most of the way there by simply measuring the chi-squared.

    Inputs
    ----------
    time
        (1d np array) Array of times.
    flux
        (1d np array) Fluxes. Median flux should be close to zero
    period, epoch, duration
        (floats) Parameters of transit. Units should be the same as the
        units of `time`. Ideally duration should not include the ingress
        and egress time. Also note the offset for `time` and `epoch` (e.g BKJD,
        TJD) should be the same.
    depth
        (float) Parameter of transit. Units should be the same as for flux
        (e.g ppm, or fractional amplitude)

    Returns
    ----------
    (float) The reduced chi squared of the fit of a box transit. A value
    of 1 indicates a perfect fit, which is implausible.

    Notes
    --------
    The expected reduced chi square value will be higher (i.e worse) for
    v-shaped transits (those that are least box shaped, and high SNR transits.
    At higher SNR, the discrepancy between the data and the box model become
    larger with respect to the typical scatter, increasing the value.
    Applying a straight cut using this metric is ill advised. Scaling the
    value by target star magnitude might work better.
    """

    idx = markTransitCadences(time, period, epoch, duration)

    #13 cadences is too long from TESS, this is a place holder
    scatter_ppm = noise.computeSgCdpp_ppm(flux[~idx], transitDuration_cadences=13)
    scatter_frac = scatter_ppm / 1e6 * np.sqrt(13)

    #Remember depth is +ve
    chi = np.sum( (flux[idx] + depth) / scatter_frac)
    chisquared = chi**2

    plt.clf()
    epoch -= 10*period
    phase = np.fmod(time - epoch + .25*period, period)
    plt.plot(phase[~idx], flux[~idx], '.', color='grey')
    plt.plot(phase[idx], flux[idx], 'k.')

    t1 = .25*period - .5*duration
    t2 = .25*period + .5*duration
    plt.plot([t1,t2], [-depth, -depth], 'r-')

    t0 = .25*period - 2.5*duration
    print(t0)
#    plt.errorbar(phase, flux, scatter_frac, color='r', lw=4, fmt='.')
    plt.errorbar(t0, 0, scatter_frac, color='r',
                 lw=2,
                 fmt='o',
                 capsize=6,
                 capthick=2,
                 zorder=+100)
#    plt.axhline(
    return chisquared / float(np.sum(idx)) - 1




