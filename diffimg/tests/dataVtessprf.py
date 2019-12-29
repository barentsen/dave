
"""
Created on Sun Dec  9 13:21:11 2018

Comparing the TESS PRF model to actual data
@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl
import pandas as pd
import numpy as np



import dave.fileio.pyfits as pyfits
import dave.diffimg.disp as disp
import dave.diffimg.tessprf
import tesscentroid


def prep():
    fn = "/home/fergal/data/tess/ffi/tess2018218052942-s0001-1-1-0120-s_ffic.fits"
    data = pyfits.getdata(fn)
    subimg = data[:848, :800]
    
    clim = np.percentile(subimg.flatten(), [2,98])
    plt.clf()
    disp.plotImage(subimg, clim=clim)

def plot1():
    fn = "/home/fergal/data/tess/ffi/tess2018218052942-s0001-1-1-0120-s_ffic.fits"
    data = pyfits.getdata(fn)
    subimg = data[:448, :400]

    c0, r0 = 56.7, 238.5
    _plot(subimg, c0, r0)
    plt.savefig('tess1-1-1-eg1.png')    

def plot2():
    fn = "/home/fergal/data/tess/ffi/tess2018218052942-s0001-1-1-0120-s_ffic.fits"
    data = pyfits.getdata(fn)
    subimg = data[:448, :400]

    c0, r0 = 184.3, 110.9
    _plot(subimg, c0, r0)
    plt.savefig('tess1-1-1-eg2.png')    


def plot3():
    fn = "/home/fergal/data/tess/ffi/tess2018218052942-s0001-1-1-0120-s_ffic.fits"
    data = pyfits.getdata(fn)
    subimg = data[:848, :800]

    c0, r0 = 564.4, 411.0
    _plot(subimg, c0, r0)
    plt.savefig('tess1-1-1-eg3.png')    




def _plot(img, c0, r0):
    prfpath = "/home/fergal/data/tess/prf"
    prfObj = dave.diffimg.tessprf.TessPrf(prfpath)
    buf = 6
    
    c1, r1 = int(c0), int(r0)
    pstamp = img[r1-buf:r1+buf+1, c1-buf:c1+buf+1]

    if False:
        prfObj = dave.diffimg.tessprf.TessPrf(prfpath)
        model = prfObj.getPrfAtColRow(c0, r0, 1, 1, 1)
        model *= np.max(pstamp) / np.max(model)
    else:
        guess = [c0, r0, np.max(pstamp), np.median(pstamp)]
        guessModel = prfObj.getPrfAtColRow(c0, r0, 1,1,1)
        guess[2] /= np.max(guessModel)


        soln = fitPrfModel(pstamp, guess, prfObj)
        print (soln)
    
        col, row, scale, sky = soln.x
        model = scale * prfObj.getPrfAtColRow(col, row, 1, 1, 1)
        model += sky
        
    plotDiagnostic(pstamp, model)
    plt.suptitle("Sector 1: Camera 1: CCD 1: Col: %g Row %g" %(c0, r0))
    

def plotDiagnostic(image, model):
    
    clim = np.log10(np.percentile(image.flatten(), (2, 98)))
    plt.clf()
    plt.subplot(131)
    plt.title("Image")
    disp.plotImage(image, log=True, clim=clim)
    
    plt.subplot(132)
    plt.title("Model")
    disp.plotImage(model, log=True, clim=clim)

    plt.subplot(133)
    diff = image - model
    plt.title("Diff Image: Score=%g" %(diffScore(diff)))
    disp.plotDiffImage(diff)


import scipy.optimize as spOpt

def fitPrfModel(image, guess, prfObj):
    """guess is col, row, amp
    """
    assert len(guess) == 4

    def callback(xk):
        pass
#        print(xk, costFunc(xk, image, prfObj)/1e6)
#        col, row, scale = xk
#        model = scale * prfObj.getPrfAtColRow(col, row, 1, 1, 1)
#        plotDiagnostic(image, model)
#        plt.pause(.1)
    
    nr, nc = image.shape
    soln = spOpt.minimize(costFunc, guess, args=(image, prfObj), 
                          method='Nelder-Mead', callback=callback)

    #These bounds are wrong because need to understand the bbox
#    bounds = [
#            (0, nc),
#            (0, nr),
#            (0, 10*np.max(image)),
#            ]            
#    soln = spOpt.minimize(costFunc, guess, args=(image, prfObj), 
#                          method='L-BFGS-B', 
#                          bounds=bounds)
    return soln


def costFunc(guess, img, prfObj):
    """This is hacky"""
    model = prfObj.getPrfAtColRow(guess[0], guess[1], 1, 1, 1)
    model *= guess[2]
    model += guess[3]
    
    diff = img - model
    cost = np.sum(diff**2)
    return cost


def diffScore(diff):
    nr, nc = diff.shape
    score = np.sum(np.fabs(diff)) / (nc*nr)
    return score