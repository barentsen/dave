# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:42:36 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import numpy as np

import dave.trapezoidFit.estimateSnr as tf
import dave.vetting.ModShift as ModShift
import dave.pipeline.clipboard as dpc
import dave.pipeline.task as task

import dave.fileio.loadMultipleDetrendings as lmd
import dave.pipeline.clipboard as clipboard
import os

def applyFunction(func, label, clip):
    """Apply function ``func`` to all detrend types"""
    detrendTypes = clip.config.detrendTypes

    clip[label] = dpc.Clipboard()
    for d in detrendTypes:
        key = '%s.%s' %(label, d)
        result = func(clip, d)
        clip[key] = result

    return clip


@task.task
def pickDefaultParams(clip):
    clip['params'] = clipboard.Clipboard()
    clip['params.period'] = 4.1591739
    clip['params.epoch'] = 2000
    clip['params.duration_hrs'] = 4
    clip['params.depth_frac'] = 400e-6
    return clip


@task.task
def serveTask(clip):
    k2id = clip['value']
    campaign = clip['config.campaign']
    storeDir = clip['config.dataStorePath']
    detrendTypes = clip['config.detrendTypes']
    assert hasattr(detrendTypes, "__len__")

    clip['serve'] = lmd.loadMultipleDetrendings(k2id, campaign, storeDir,
                         detrendTypes)


    #Enforce contract. (Make sure expected keys are in place)
    clip['serve.time']
    clip['serve.cube']
    clip['serve.%s.flux' %(detrendTypes[0])]
    clip['serve.tpfHeader']
    return clip


@task.task
def trapFitTask(clip):
    clip = applyFunction(trapFitWrapper, 'trapFit', clip)
    return clip


def trapFitWrapper(clip, detrendKey):

    time_days = clip['serve.time']
    flux_norm = clip['serve.%s.flux' %(detrendKey)]
    flags = clip['serve.%s.flags' %(detrendKey)]
    period_days = clip['params.period']
    duration_hrs = clip['params.duration_hrs']
    phase_bkjd = clip['params.epoch']  #Check this what BLS returns
    depth_frac = clip['params.depth_frac']


    #We don't know these values.
    unc = np.ones_like(flux_norm)
    unc[flags] = 1e99
    flux_norm[flags] = 0

    assert(np.all(np.isfinite(time_days[~flags])))
    assert(np.all(np.isfinite(flux_norm[~flags])))

    try:
        out = tf.getSnrOfTransit(time_days, flux_norm,\
            unc, flags, \
            period_days, phase_bkjd, duration_hrs, depth_frac)
    except ValueError, e:
        print "trap fit failed. Stop ignoring this errror"
        return clip

    #Enforce contract
    out['period_days']
    out['epoch_bkjd']
    out['duration_hrs']
    out['ingress_hrs']
    out['depth_frac']
    out['bestFitModel']
    out['snr']
    return out


from dave.plot.compare_lightcurves import plot_multiple_lightcurves
@task.task
def makeOverviewPlotsTask(clip):

    detrendTypes = clip.config.detrendTypes

    epic = clip.value
    campaign = clip.config.campaign
    period = clip.params.period
    epoch = clip.params.epoch
    basePath = clip.config.onepageBasename

    time = clip.serve.time
    fluxList = []

    for d in detrendTypes:
        flux = clip['serve.%s.flux' %(d)]
        flags = clip['serve.%s.flags' %(d)]

        fluxList.append(flux)

    fluxList = tuple(fluxList)
    fileRoot = "%i-%i-%i.png" %(epic, campaign, int(period*10))

#    mp.clf()
#    for f in fluxList:
#        mp.plot(time, f, 'o', ms=2)
    title="EPIC %i Campaign %i" %(epic, campaign)
    fig = plot_multiple_lightcurves(time, fluxList, labels=detrendTypes,
                              title=title, mec="none")
#    fn = os.path.join(basePath, "lc"+fileRoot)
#    fig.savefig(fn)
#
    title="EPIC %i Campaign %i: Period: %.4f (days) Epoch: %.3f" \
            %(epic, campaign, period, epoch)
    plot_multiple_lightcurves(time, fluxList, labels=detrendTypes,
                              title=title, period=period, epoch=epoch)
#    fn = os.path.join(basePath, "fold"+fileRoot)
#    fig.savefig(fn)

    return clip


