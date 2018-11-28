#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Create tasks to run a vetting pipeline on tess data.
Created on Tue Nov 27 21:23:14 2018

@author: smullally

"""

import numpy as np

 

@task.task
def lppMetricTask(clip):
    
    class clipToLppInputClass(object):
    
        def __init__(self, clip):
            """
            create a TCE class from the clipboard info
            """
            self.time=clip['serve.time']
            self.tzero=clip['serve.param.epoch_btjd']
            self.dur=clip['serve.param.transitDuration_hrs']
            self.period=clip['serve.param.orbitalPeriod_days']
            self.flux=clip['serve.detrendFlux']
                    
            
    data=clipToLppInputClass(clip)
    mapInfo=clip['config.lppMapFilePath']
    normTLpp,rawTLpp,transformedTransit=computeLPPTransitMetric(data,mapInfo)

    out = dict()
    out['TLpp'] = normTLpp
    out['TLpp_raw'] = rawTLpp

    clip['lpp'] = out

    #Enforce contract
    clip['lpp.TLpp']

    return clip


import dave.vetting.ModShift as ModShift
@task.task
def modshiftTask(clip):
    
    time = clip['serve.time']
    flux = clip['detrend.detrendFlux']
    fl = np.zeros(len(time),dtype=bool)
    period_days = clip['trapFit.orbitalPeriod_days']
    epoch_btjd = clip['trapFit.epoch_btjd']
    detrendType = clip.config.detrendType

    tic = clip['config.tic']
    basePath = clip['config.modshiftBasename']
    ticStr = "%016i" %(tic)
    basename = getOutputBasename(basePath, ticStr)

    # Name that will go in title of modshift plot
    objectname = "TIC %012i" % (tic)

    model = clip['serve.modelFlux']
    modplotint=int(1)  # Change to 0 or anything besides 1 to not have modshift produce plot
    plotname = "%s-%02i-%04i-%s" % (basename,np.round(clip.bls.period*10),np.round(clip.bls.epoch),clip.config.detrendType)
    out = ModShift.runModShift(time[~fl], flux[~fl], model[~fl], \
        plotname, objectname, period_days, epoch_btjd, modplotint)
    clip['modshift'] = out

    #Enforce contract
    clip['modshift.mod_Fred']
    clip['modshift.mod_ph_pri']
    clip['modshift.mod_secdepth']
    clip['modshift.mod_sig_pri']
    return clip