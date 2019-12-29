# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 06:39:07 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

#211923431.PC
#depth: 0.000854962
#duration_hrs: 4
#epoch: 2310.83
#period: 29.7152
#snr: 4.8528
#BLS SNR= 6.21
#
#
#211945775.0  weak PC
#depth: 0.00160972
#duration_hrs: 6
#epoch: 2309.63
#period: 25.0146
#snr: 2.89823
#BLS SNR= 3.02
#
#1 211911246.0 weak PC
#bls_search_periods: np.ndarray (7254,)
#convolved_bls: np.ndarray (7254,)
#depth: 0.000835716
#duration_hrs: 6
#epoch: 2309.35
#period: 26.7379
#snr: 2.94035


def exporterTask(clip):
    """Export a text file of most important properties

    TODO:
    Should output to path suggested by clipboard
    """
    epic = clip['value']
    campaign = clip['config.campaign']
    fit = clip['trapFit']
    disp = clip['disposition']

    hdr = ["Epic"]
    vals = ["%9i" % epic]

    keys = "period_days epoch_bkjd duration_hrs depth_frac snr".split()
    hdr.extend(keys)
    for k in keys:
        vals.append("%10.3f" %(fit[k]))

    hdr.append("RpRs")
    vals.append( "%10.3e" %(np.sqrt(fit['depth_frac'])))

    keys = "isCandidate isSignificantEvent".split()
    hdr.extend(keys)
    for k in keys:
        vals.append("%i" %(disp[k]))

    hdr.append("reasonForFail")
    vals.append(disp['reasonForFail'])

    out = ["#" + " |".join(hdr)]
    out.append(" |".join(vals))

    fn = "k%09i-c%02i.export" %(epic, campaign)
    fp = open(fn, 'w')
    fp.write("\n".join(out))
    fp.close()