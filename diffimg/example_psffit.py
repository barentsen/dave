"""
Created on Wed Nov 14 21:34:47 2018

Example code for psffit.py
@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import numpy as np

import dave.fileio.mastio as kmastio
import dave.fileio.tpf as ktpf
import datetime
import psffit


def timing(numIter):
    """Function to help measure timings"""

    guess = [3., 3.5, 1.1, 1000, 10]
    nc, nr, = 9, 12

    func = psffit.gaussianWithConstantSkyPsf
#    func = wrapper
    t0 = datetime.datetime.now()

    for i in xrange(numIter):
        guess = [3.+ .01*i, 3.5- .01*i, 1.1, 1000, 10]
        psffit.computeModel(nc , nr , func, guess)

    t1 = datetime.datetime.now() - t0
    print(t1)
    sec = t1.total_seconds()
    print ("Took %.3f sec, %.3f msec/itr" %(sec, 1e3*sec/float(numIter)))





def fitKepler():
    """Fit some Kepler data as an example"""
    kepid = 8311864
    quarter = 6

    ar= kmastio.KeplerArchive()
    fits, hdr = ar.getLongTpf(kepid, quarter, header=True)
    cube = ktpf.getTargetPixelArrayFromFits(fits, hdr)

    out = []
    for i in range(200, 201):
        print("Image %i" %(i))
        img = cube[i,:,:]
        nr, nc = img.shape

        func = psffit.gaussianWithConstantSkyPsf
        guess = (3.0, 1.8, .41, 3e4, -2000)
        soln = psffit.fitPrf(img, func, guess)
        print (soln)
        out.append(soln.x)

        plt.clf()
        model = psffit.computeModel(nc, nr, func, soln.x)

        plt.subplot(131)
        plt.title("Data")
        plt.imshow(np.log10(img), origin='bottom')
        plt.colorbar()

        plt.subplot(132)
        plt.title("Model")
        plt.imshow(model, origin='bottom')
        plt.colorbar()

        plt.subplot(133)
        plt.title("Diff")
        plt.imshow(img - model, origin='bottom')
        plt.colorbar()

        plt.suptitle("%i Q%i RIN %i" %(kepid, quarter, i))
        plt.savefig('k%i-Q%i-R%i.png' %(kepid, quarter, i))
        plt.pause(.1)

    out = np.array(out)


if __name__ == "__main__":
    timing(2)