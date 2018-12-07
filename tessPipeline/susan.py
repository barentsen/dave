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
import dave.tessPipeline.tessPipeline as tessPipeline
import dave.vetting.RoboVet as RoboVet
import pandas as p

def createConfig(sector, tic, planetNum, debugMode=True):
    cfg = dict()
    cfg['debug'] = debugMode
    cfg['sector'] = sector
    cfg['tic'] = tic
    cfg['planetNum'] = planetNum
    cfg['value'] = tic

    #TODO This shouldn't be hardcoded, but passed as a parameter
    cfg['dvtLocalPath'] = "/Users/smullally/TESS/TCEs/Sector1/dvt/"
    
    #TODO Need modshift paths
    cfg['lppMapFile'] = "/Users/smullally/Code/lpptransitlikemetric/data/maps/combMapDR25AugustMapDV_6574.mat"

    cfg['modshiftBasename'] = "/Users/smullally/TESS/TCEs/Sector1/dave"    
    
    cfg['taskList'] = ['serveTask', 'lppMetricTask', 'modshiftTask','sweetTask']
    
    clip = clipboard.Clipboard(cfg)
    
    return clip
    

def outputInfo(clip):
    """Createa  text string to output to a file
    """     
    header="tic,planetNum"
    
    text = "%i" % (clip.config.tic)
    text = "%s,%i" % (text, clip.config.planetNum)
    
    for k in clip.serve.param.keys():
        text = "%s,%f" % (text, clip['serve']['param'][k])
        header = "%s,%s" % (header, k)

    text = "%s,%f" % (text, clip.lpp.TLpp)
    header = "%s,Tlpp" % (header) 
    
    for k in clip.modshift.keys():    
        text = "%s,%f" % (text, clip['modshift'][k])
        header = "%s,%s" % (header, k)
    
    for k in clip.robovet.keys():
        text = "%s, %s" % (text, str(clip['robovet'][k]))
        header = "%s,rv_%s" % (header, k)
    
    return text,header


def runOneDv(sector,tic,planetNum,debugMode=True):
    
    cfg = createConfig(sector,tic,planetNum,debugMode=True)
    
    clip = tessPipeline.runOne(cfg,returnClip=True)
    
    # This is an attempt at quickly vetting the signals.
    # This should be its own task.
    
    rv = RoboVet.roboVet(clip.modshift)
    
    lpp_th = 4  #Threshold for LPP
    
    if clip.lpp.TLpp > lpp_th:
        rv['disp'] = 'false positive'
        rv['comments']= rv['comments'] + "-LPP_TOO_HIGH"
        rv['not_trans_like']=1
    
    sweet_th = 3.5
    if (clip.sweet.amp[0,-1] > sweet_th) | \
       (clip.sweet.amp[1,-1] >sweet_th) | \
       (clip.sweet.amp[2,-1] > sweet_th):
           rv['disp']='false positive'
           rv['comments']=rv['comments'] + "-SWEET_FAIL"
           rv['not_trans_like']=1
    
    clip['robovet']=rv
    

    return clip
    
    
    
def runAllTces(tceFile,sector,outfile):
    """
    Run for all TCEs
    """
    df=p.read_csv(tceFile,comment='#')
    
    for index,row in df[0:5].iterrows():
        
        clip=runOneDv(sector,row.ticid,row.planetNumber)
        text,header = outputInfo(clip)
        
        #Write out decision
        with open(outfile,"a") as fp:
        
            if index == 0:
                fp.write(header+"\n")
                
            fp.write(text + "\n")
    return clip
    