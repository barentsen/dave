# -*- coding: utf-8 -*-
"""
TODO, use requests to download data from http as needed.
"""

from pdb import set_trace as debug

from abstracttessprf import AbstractTess

import astropy.io.fits as pyfits
import numpy as np
import requests
import os

class TessFitsPrf(AbstractTess):

    def __init__(self, cache_path, sector, camera, ccd):
        AbstractTess.__init__(self, cache_path, sector, camera, ccd)

        #This location must exist, otherwise who knows where we will
        #end up storing the downloaded fits files
        assert os.path.exists(cache_path)

        #Once we load an image from disk, keep it in memory for faster
        #access
        self.imgCache = dict()

        self.url = 'https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/'
        #How many model evalulations are there in the column direction?
        self.modelsPerCol = 5

        #Model locs encodes the *column* and *row* (in that order)
        #of the locations where the PRF models are realised
        #The class will assume this variable encodes a regular grid
        #of models in col and row
        self.modelLocs = np.array([ \
                                       [45, 1], \
                                       [557, 1], \
                                       [1069, 1], \
                                       [1580, 1], \
                                       [2092, 1], \
                                       [45, 513], \
                                       [557, 513], \
                                       [1069, 513], \
                                       [1580, 513], \
                                       [2092, 513], \
                                       [45, 1025], \
                                       [557, 1025], \
                                       [1069, 1025], \
                                       [1580, 1025], \
                                       [2092, 1025], \
                                       [45, 1536], \
                                       [557, 1536], \
                                       [1069, 1536], \
                                       [1580, 1536], \
                                       [2092, 1536], \
                                       [45, 2048], \
                                       [557, 2048], \
                                       [1069, 2048], \
                                       [1580, 2048], \
                                       [2092, 2048], \
                                   ])

        #Check previous two variables are consistent
        assert len(self.modelLocs) % self.modelsPerCol == 0
        #Check Row numbers are always ascending
        assert np.all(np.diff(self.modelLocs[:,1]) >= 0)


    def getPrfAtColRow(self, col, row, *args):
        assert len(args) == 0
        col = float(col)
        row = float(row)

        brCol, brRow = self.getBracketingPixelPositions(col, row)
        subSampledImgList = self.getBracketingImages(brCol, brRow)

        subSampledModel = self.interpolatePrf(subSampledImgList,
                                              col, row,
                                              brCol, brRow)


        regSampledModel =  self.getRegularPrfFromSubsampledPrf(subSampledModel,
                                                               col, row)

        return regSampledModel


    def getBracketingPixelPositions(self, col, row):

        """

        modelLocs = [[c0, r0], [c1, r0], ...]

        Returns
        ---------
        brCol
            (np array) Array of column positions of the 4 model realisations
            that bracket the requested col, row

        brRow
            (np array)  Array of row positions of the 4 model realisations
            that bracket the requested col, row


        brCol is in the order [c0, c0, c1, c1]
        brRow is in the order [r0, r1, r0, r1], according to this diagram,
        (with c0 < c1, and r0 < r1)

             (c0, r1)    (c1, r1)
            +-----------+
            |           |
            |    (c, r) |
            |   o       |
            |(c0,r0)    | (c1, r0)
            +-----------+

        """

        modelColLocations = self.modelLocs[:,0]
        modelRowLocations = self.modelLocs[:,1]

        #Note, these assumptions are valid for TESS, but not for a more
        #general method (e.g for Kepler) some valid locations defy these
        #constraints
        assert col >= modelColLocations[0]
        assert col < modelColLocations[-1]
        assert row >= modelRowLocations[0]
        assert row < modelRowLocations[-1]

        idx = (col >= modelColLocations) & (row >= modelRowLocations)
        wh = np.where(idx)[0][0]
        c0, r0 = self.modelLocs[wh]
        c1, r1 = self.modelLocs[wh + self.modelsPerCol + 1]
        assert c1 > c0
        assert r1 > r0

        brCol = np.array([c0, c1, c0, c1])
        brRow = np.array([r0, r0, r1, r1])

        return brCol, brRow


    def getBracketingImages(self, brCol, brRow):
        fitsFileList = self.getFitsFilenames(brCol, brRow)
        subSampledImgList = self.loadFitsFiles(fitsFileList)
        return subSampledImgList



    def getFitsFilenames(self, brCol, brRow):
        """Deduce fits file names for the input bracketing columns and rows

        Makes no effort to ensure the filenames are valid or that the files
        exist
        """

        assert len(brCol) == len(brRow)
        #This is the portion of the path (excluding the filename) that
        #is shared between the remote url and the cached location on disk
        imageSubPath = os.path.join("start_s%04i" %(self.sector),
                            "cam%i_ccd%i" %(self.camera, self.ccd))

        fitsDateStr = self.getFitsDateStr(self.sector)

        f = lambda x, y:  "tess%s-prf-%i-%i-row%04i-col%04i.fits" \
                            %(fitsDateStr,
                              self.camera,
                              self.ccd,
                              y, x)  #names are -row-col.fits

        flist = list(map(f, brCol, brRow))
        flist = list(map(lambda x: os.path.join(imageSubPath, x), flist))
        return flist


    def getFitsDateStr(self, sector):
        """Mapping from sector number to fits filename's datestring

        PRF model fits files include a date string whose meaning is unclear,
        and uncomputable. We need a look up table, which is stored here.

        Models are only generated for select sectors. Use sectorLookup()
        in the base Tess class to convert your actual sector to the sector
        number of the PRF files you need.
        """
        if sector == 4:
            return '2019107181900'
        return '2018243163600'


    def loadFitsFiles(self, flist):
        out = []

        for f in flist:
            out.append(self.loadSingleImage(f))
        return out


    def loadSingleImage(self, imageSubPath):
        """
        Look for the image first in memory, then on disk, then on the web.

        The image should be expected to exist at both
        self.path + imageSubPath and self.url + imageSubPath

        Returns
        ------
        A numpy 2d array
        """

        key = imageSubPath
        if key not in self.imgCache:
            cache_path = os.path.join(self.path, imageSubPath)

            if not os.path.exists(cache_path):
                url = os.path.join(self.url, imageSubPath)
                self.download(url, cache_path)

            self.imgCache[key] = pyfits.getdata(cache_path)

        img = self.imgCache[key]
        assert img.ndim == 2
        return img


    def download(self, remoteUrl, local):
        """Download the file at `remoteUrl` to the path `local` on disk"""

        localDir = os.path.split(local)[0]
        if not os.path.exists(localDir):
            os.makedirs(localDir)

        print(remoteUrl)
        debug()
        r = requests.get(remoteUrl, allow_redirects=True)
        open(local, 'wb').write(r.content)
