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
    flux = clip['detrend.flux_frac']
    fl = clip['detrend.flags']
    period_days = clip['trapFit.period_days']
    epoch_bkjd = clip['trapFit.epoch_bkjd']
    detrendType = clip.config.detrendType
    #dur_hrs =  clip['trapFit.duration_hrs']
    #ingress_hrs = clip['trapFit.ingress_hrs']
    #depth_ppm = 1e6*clip['trapFit.depth_frac']

    epic = clip['value']
    basePath = clip['config.modshiftBasename']
    epicStr = "%09i" %(epic)
    basename = getOutputBasename(basePath, epicStr)

    # Name that will go in title of modshift plot
    objectname = "EPIC %09i" %(epic)
#
#    subSampleN= 15
#    ioBlock = trapFit.trapezoid_model_onemodel(time[~fl], period_days, \
#        epoch_bkjd, depth_ppm, dur_hrs, \
#        ingress_hrs, subSampleN)
#    model = ioBlock.modellc -1   #Want mean of zero
#    #model *= -1  #Invert for testing
    model = clip['trapFit.bestFitModel']
    modplotint=int(1)  # Change to 0 or anything besides 1 to not have modshift produce plot
    plotname = "%s-%02i-%04i-%s" % (basename,np.round(clip.bls.period*10),np.round(clip.bls.epoch),clip.config.detrendType)
    out = ModShift.runModShift(time[~fl], flux[~fl], model[~fl], \
        plotname, objectname, period_days, epoch_bkjd, modplotint)
    clip['modshift'] = out

    #Enforce contract
    clip['modshift.mod_Fred']
    clip['modshift.mod_ph_pri']
    clip['modshift.mod_secdepth']
    clip['modshift.mod_sig_pri']
    return clip