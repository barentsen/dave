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



from dave.task import task

@task
def serveTask(clip):

    pass

    #Enforce contract
    clip['serve.time']
    clip['serve.detrend_flux']
    clip['serve.param.orbital_period_days']
    clip['serve.param.epoch_btjd']
    clip['serve.param.transitDepth_ppm']
    clip['serve.param.transitDuration_hrs']     

    return clip