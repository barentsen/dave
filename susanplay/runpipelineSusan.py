# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:20:23 2015

@author: sthomp
"""

import dave.susanplay.mainSusan as mS
#import dave.pipeline.pipeline as pipe
import numpy as np
#import dave.pipeline.plotting as pp
import matplotlib.pyplot as plt
import dave.susanplay.sueplotting as sp
import dave.plot.multipage as mp

    
#%%   
#infile='/home/sthomp/DAVE/playK2/k2_go3049.txt'
vers="20150128";
infile='/home/sthomp/DAVE/playK2/KEES_2016_01_13.txt'
outfile='/home/sthomp/DAVE/playK2/KEES_2016_01_13_out%s.txt' %(vers)
#outcand='/home/sthomp/DAVE/playK2/k2_list_cand.txt'
fid=open(outfile,'a')

cfg = mS.loadMyConfiguration()
cfg['debug'] = False
cfg['modshiftBasename']='/home/sthomp/DAVE/playK2/Kees/k';
#cfg['prfPath']='morejunk/junk';

indata=np.loadtxt(infile,dtype='float',delimiter=None,comments='#',usecols=[0,5])
want=(indata[:,1]>=3) & (indata[:,1] <=5)
data=indata[want,:]
#%%
span=[55,100]
for i,v in enumerate(data[span[0]:span[1],0]):    
    epicid=np.int(v)
    acti=i+span[0]
    cfg['campaign'] = int(data[acti,1])
    print epicid,cfg['campaign']

    output=mS.runOne(epicid, cfg)
     
    try:
        print output['exception']
        rep="%u --bad -- %s\n" % (epicid, output['exception'])
        fid.write(rep)
        print ' ! !  EXCEPTION EXCEPTION ! ! ! ! '
        oneout="/home/sthomp/DAVE/playK2/epic%u.exc" % epicid
        fida=open(oneout,'w')
        fida.write(output['backtrace'])
        fida.close()
        
    except KeyError:
        
        try:
        
#            plt.figure(1)
#            sp.summaryPlot(output)
#            outfig="%sfig%s.png" % (cfg['modshiftBasename'],str(epicid))
#            plt.savefig(outfig)
#            plt.figure(2)
#            sp.indivPlot(output,6)
#            outfig="%sind%s.png" % (cfg['modshiftBasename'],str(epicid))
#            plt.savefig(outfig)
#            plt.pause(.1)
            if output.disposition.isCandidate == 1:
                disp="CANDIDATE"
            else:
                disp="FALSE POSITIVE"
            output.disposition.finaldisp=disp
            
            period=output.bls.period
            #info = "%u  %6.2f  %s   %s\n" % (epicid, period, disp, output['disposition.reasonForFail'])
            line, hdr= createExportString(output,delimiter="|")        
            fid.write("%s  |%s\n" % (line, disp))
            
            outfile="%s%s-%smp.pdf" % (cfg['modshiftBasename'],str(epicid),vers)
            info="version %s\n%u\n%s\n%s" % (vers,epicid,disp,output.disposition.reasonForFail)
            mp.plot_all_multipages(outfile,output,info)
            print outfile
            
        except KeyError, IOError:
            print epic
            print "Could not print outfile"
        
fid.close()

#%%
#For running known candidates
vers='20160205'
ephemfile='k2candidatesc3.csv';
outfile='k2candidatesc3_out.txt';
fid=open(outfile,'a')

cfg = mS.loadMyConfiguration()
cfg['debug'] = False
cfg['modshiftBasename']='/home/sthomp/DAVE/playK2/pcs/';
#cfg['prfPath']='morejunk/junk';


indata=np.loadtxt(ephemfile,dtype='string',delimiter=',',comments='#',usecols=[2,5,6,10,14,15])

periods=indata[:,2].astype(float)
ids=np.floor(indata[:,0].astype(float))
epochs=indata[:,3].astype(float)-2454833;
depths=indata[:,4].astype(float)/100.0;
durs=indata[:,5].astype(float)*24.0;
camps=indata[:,1].astype(int)
#%
span=[0,len(periods)]
for i,v in enumerate(ids[span[0]:span[1]]):    
    epicid=int(v)
    acti=i+span[0]
    cfg['campaign'] = int(camps[acti])
    print epicid,cfg['campaign'], periods[acti]

    output=mS.runOneEphem(ids[acti],periods[acti],epochs[acti], cfg,durs[acti],depths[acti])


    try:
        
        print output['exception']
        fid.write(('%f  %f  Exception\n' % (epicid,periods[acti])))
        oneout="/home/sthomp/DAVE/playK2/pcs/epic%u%u.exc" % (epicid,np.floor(periods[acti]))
        fida=open(oneout,'w')
        fida.write(str(output['exception']))
        fida.write(output['backtrace'])
        fida.close()
        
    except KeyError:        
        if output.disposition.isCandidate == 1:
            disp="CANDIDATE"
        else:
            disp="FALSE POSITIVE"
            
        output.disposition.finaldisp=disp    
        line, hdr= mp.createExportString(output,delimiter="|")        
        fid.write("%s  |%s|  %f\n" % (line, disp, epochs[acti]))
        plt.figure(10)
        sp.summaryPlot(output)
        try:
            outfile="%s%s-%s-%smp.pdf" % (cfg['modshiftBasename'],str(epicid),vers,np.floor(periods[acti]))
            info="version %s\n%u\n%s\n%s" % (vers,epicid,disp,output.disposition.reasonForFail)
            mp.plot_all_multipages(outfile,output,info)
        except:
            pass
    
fid.close()
