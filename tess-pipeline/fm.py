# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Tue Nov 27 21:25:04 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import numpy as np


import dave.tess_pipeline.tessfunc as tessfunc
import dave.pipeline.clipboard as clipboard
from dave.pipeline.task import task


@task
def serveTask(clip):
    
    sector = clip['config.sector']
    tic = clip['config.tic']
    planNum = clip['config.planetNum']
    localPath = clip['config.dvtLocalPath']
    
    dvt, hdr = tessfunc.serve(sector, tic, planNum, localPath)

    out = dict()
    out['time'] = dvt['TIME']
    out['detrend_flux'] = dvt['DETREND_FLUX']
    
    out['param.orbitalPeriod_days'] = hdr['TPERIOD']
    out['param.epoch_btjd'] = hdr['TEPOCH']
    out['param.transitDepth_ppm'] = hdr['TDEPTH']
    out['param.transitDuration_hrs'] = hdr['TDUR']
    
    clip['serve'] = out
 
    #Enforce contract
    clip['serve.time']
    clip['serve.detrend_flux']
    clip['serve.param.orbitalPeriod_days']
    clip['serve.param.epoch_btjd']
    clip['serve.param.transitDepth_ppm']
    clip['serve.param.transitDuration_hrs']     

    return clip



