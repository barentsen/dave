# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:39:30 2016

@author: fergal

$Id$
$URL$
"""

__version__ = "$Id$"
__URL__ = "$URL$"



import matplotlib.pyplot as mp
import matplotlib as mpl
import numpy as np

import dave.fileio.kplrfits as kf
import dave.fileio.mastio as mastio2
import dave.fileio.tpf as tpf
import diffimg
import plotTpf
import arclen

def main():
    ar = mastio2.K2Archive()
    kepid = 206103150

    fits = ar.getLongCadence(kepid, 3)
    time = fits['TIME']
    flags = fits['SAP_QUALITY']
    pa = fits['SAP_FLUX']
    pdc = fits['PDCSAP_FLUX']
    cent1 = fits['MOM_CENTR1']
    cent2 = fits['MOM_CENTR2']
    badIdx = np.isnan(time) | np.isnan(pdc)

    pa /= np.nanmedian(pa)
    pa -= 1
    pa *= 1e3

    pdc /= np.nanmedian(pdc)
    pdc -= 1
    pdc *= 1e3

    #Compute roll phase
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = arclen.computeArcLength(centColRow, badIdx)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    #Load TPF
    fits, hdr = ar.getLongTpf(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)
    gain = hdr['gain']
    cube *= gain

    #Cut to just a short sectoin of interest
    idx = (time > 2207.2) & (time < 2208.8)
#    idx = np.ones_like(time, dtype=bool)
    time = time[idx]
    pa = pa[idx]
    pdc = pdc[idx]
    cube = cube[idx]
    flags = flags[idx]
    rollPhase = rollPhase[idx]


    i0 = 53
    mp.style.use('apj')
    rc = mpl.rcParams
    rc['ytick.labelsize'] = 12
    rc['axes.labelsize'] = 20
#    diffimg.plotRollPhaseDiagnosticPlot(None, rollPhase, flags, i0)

    drawFrame(cube, hdr, time, pa, pdc, rollPhase, flags, i0)
#    return


    #17
    for i in range(17, len(time) - 16):
        drawFrame(cube, hdr, time, pa, pdc, rollPhase, flags, i)
        mp.pause(.01)
        mp.savefig("wasp69-%03i.png" %(i))

def drawFrame(cube, hdr, time, pa, pdc, rollPhase, flags, i0):
    mp.clf()

    assert(len(time) == len(cube))
    assert(len(time) == len(pdc))
    assert(len(time) == len(rollPhase))
    assert(len(time) == len(flags))

#    index = np.arange(len(cube))
    gs = (12,3)
    mp.subplot2grid(gs, (0,0), rowspan=2, colspan=3)
    mp.plot(time[i0], pdc[i0], 'wo', ms=12, mec="r", mew=2)
    mp.plot(time[i0], pa[i0], 'wo', ms=12, mec="r", mew=2)
    mp.plot(time, pa, 'bo-')
    mp.plot(time, pdc, 'ko-')
    mp.ylabel("Flux (ppk)")

    mp.axvspan(2208.14, 2208.28, color='#AAAAFF', alpha=.4)
    #Draw a vertical line
    y1, y2 = mp.ylim()
    y3 = max(pa[i0], pdc[i0])
    mp.plot([time[i0], time[i0]], [y1, y3], color='grey')

    mp.subplot2grid(gs, (3,0), rowspan=2, colspan=3)
    mp.plot(time[i0], rollPhase[i0], 'wo', ms=12, mec="r", mew=2)
    mp.axhline(rollPhase[i0], color='r', lw=1)
    mp.plot(time, rollPhase, 'ko')
    mp.ylabel("Roll angle")
    mp.ylim(-1,1)
    mp.axvspan(2208.14, 2208.28, color='#AAAAFF', alpha=.4)

    #Draw a vertical line
    mp.plot([time[i0], time[i0]], [rollPhase[i0], 1], color='grey')

    diff, oot, diagnostics = \
        diffimg.constructK2DifferenceImage(cube,  i0, rollPhase, flags)

    if diagnostics['errorMsg'] != "None":
        diff *= 0

    print i0, diagnostics
    mp.subplot2grid(gs, (6,0), rowspan=6, colspan=1)
    plotTpf.plotCadence(cube[i0], hdr)
    mp.clim(0, 1.20e7)
    mp.colorbar()
    mp.xlabel("Cadence Image")
    mp.gca().set_xticks([])
    mp.gca().set_yticks([])

    mp.subplot2grid(gs, (6,1), rowspan=6, colspan=1)
    plotTpf.plotCadence(diff, hdr)
    mp.clim(-1e4, 1e4)
    mp.colorbar()
    mp.xlabel("Difference Image")
    mp.gca().set_xticks([])
    mp.gca().set_yticks([])
#
    mp.subplot2grid(gs, (6,2), rowspan=6, colspan=1)
    snr  = diff / np.sqrt(np.fabs(cube[i0]))
    plotTpf.plotCadence(snr, hdr)
    mp.colorbar()
    mp.clim(0, 30)
    mp.xlabel("Diff. SNR")
    mp.gca().set_xticks([])
    mp.gca().set_yticks([])


def plot2():
    mp.style.use('apj')
    ar = mastio2.K2Archive()
    kepid = 206103150

    fits = ar.getLongCadence(kepid, 3)
    time = fits['TIME']
    flags = fits['SAP_QUALITY']
    pa = fits['SAP_FLUX']
    pdc = fits['PDCSAP_FLUX']
    cent1 = fits['MOM_CENTR1']
    cent2 = fits['MOM_CENTR2']
    badIdx = np.isnan(time) | np.isnan(pa)
    time[ badIdx ] = 0

#    return flags

    #Compute roll phase
    centColRow = np.vstack((cent1, cent2)).transpose()
    rot = arclen.computeArcLength(centColRow, badIdx)
    rollPhase = rot[:,0]
    rollPhase[flags>0] = -9999    #A bad value

    mp.clf()
    x = np.arange(len(rollPhase))
    mp.plot(time, rollPhase, 'ko', label="K2 Roll angle", ms=14)

    #Plot thruster firings
    thrusterFiring = (flags & kf.SapQuality['DefiniteRollTweak']) > 0
    idx2 = np.roll(thrusterFiring, -1)
    mp.plot(time[thrusterFiring], rollPhase[idx2], 'r*', ms=36, label="Thruster Firing")


    #Mark a single cadence
    i0 = np.argmin( np.fabs(time - 2168.39) )
    y = rollPhase[i0]
    mp.plot([2168.15, 2168.73], [y,y], 'g-')
    mp.plot(time[i0], rollPhase[i0], 'gs', ms=20, mew=0)

    #Mark cadences to be interpolated
    s = [i0-13, i0-12, i0+16, i0+17]
    mp.plot(time[s], rollPhase[s], 'mo', label="Reference Cadences", ms=14)


#
    mp.axis([2167.457, 2169.46, -.54, +.9])
    mp.xlabel("Time (BKJD)")
    mp.ylabel("Roll Angle")
    mp.legend(loc=4)
    ax = mp.gca()
    for a in [ax.xaxis, ax.yaxis]:
        if a.get_scale() == 'linear':
            a.set_minor_locator(mpl.ticker.AutoMinorLocator())

        a.set_tick_params(which="major", length=8, width=1, color='k')
        a.set_tick_params(which="minor", length=4, width=1, color='k')



    mp.savefig('example-roll1.png')
    if True:
        #Mark a single cadence that I can't make difference images
        i0 = np.argmin( np.fabs(time - 2168.62) )
        y = rollPhase[i0]
        mp.plot([2168.31, 2168.88], [y,y], 'g-')
        mp.plot(time[i0], rollPhase[i0], 'gs', ms=20, mew=0)
        mp.text(2168.24, y-.03, "?", color='g', fontsize=48)

        s = [i0+13, i0+14]
        mp.plot(time[s], rollPhase[s], 'mo', ms=14)

    mp.savefig('example-roll2.png')



def exampleImages():
    ar = mastio2.K2Archive()
    kepid = 206103150

    #Load TPF
    fits, hdr = ar.getLongTpf(kepid, 3, header=True)
    cube = tpf.getTargetPixelArrayFromFits(fits, hdr)

    mp.clf()
    for i in [1000, 1001]:
        plotTpf.plotCadence(cube[i], hdr)
        mp.clim(0, 15e3)
        mp.gca().set_xticks([])
        mp.gca().set_yticks([])
        mp.savefig('exampleTpfImg-%i.png' %(i))