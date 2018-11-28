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


import dave.tess_pipeline.tessfunc as tessfunc
import dave.pipeline.clipboard as clipboard
from dave.pipeline.task import task


def createConfig(sector, tic, planetNum):
    cfg = dict()
    cfg['sector'] = sector
    cfg['tic'] = tic
    cfg['planetNum'] = planetNum
    
    cfg['taskList'] = ['serveTask']
    
    clip = clipboard.Clipboard(cfg)
    return clip
    
        

def runOne(config, returnClip=False):
    """Run the pipeline on a single target.

    Inputs:
    ------------
    k2id
        (int) Epic id of target to run on.

    config
        (dict) Dictionary of configuration parameters as created by, e.g
        loadMyConfiguration()

    Returns:
    ---------
    A clipboard containing the results.

    Notes:
    ---------
    Don't edit this function. The pipeline can recover gracefully from
    errors in any individual task, but an error in this function will crash
    the pipeline
    """

    taskList = config['taskList']

    clip = clipboard.Clipboard()
    clip['config'] = config

    #Check that all the tasks are properly defined
    print( "Checking tasks exist")
    for t in taskList:
        f = eval(t)

    #Now run them.
    for t in taskList:
        print( "Running %s" %(t))
        f = eval(t)
        clip = f(clip)

    if returnClip:
        return clip


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



