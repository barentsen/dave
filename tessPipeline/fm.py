# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Tue Nov 27 21:25:04 2018

astropy
scikit-learn
lpproj (pip)

@author: fergal
"""

#from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug


import dave.tessPipeline.tessfunc as tessfunc
import dave.pipeline.clipboard as clipboard
from dave.pipeline.task import task


def createConfig(sector, tic, planetNum, debugMode=True):
    cfg = dict()
    cfg['debug'] = debugMode
    cfg['sector'] = sector
    cfg['tic'] = tic
    cfg['planetNum'] = planetNum

    #TODO This shouldn't be hardcoded, but passed as a parameter
    cfg['dvtLocalPath'] = "/home/fergal/data/tess/dvt/sector1"
    
    #TODO Need modshift paths
    cfg['lppMapFile'] = "/home/fergal/data/tess/lppmaps/combMapDR25AugustMapDV_6574.mat"

    cfg['modshiftBasename'] = "/home/fergal/data/tess/daveOutput/"    
    
#    cfg['taskList'] = ['serveTask', 'lppMetricTask', 'modshiftTask']
    cfg['taskList'] = ['serveTask', 'sweetTask']
    
    clip = clipboard.Clipboard(cfg)
    return clip
    
        
