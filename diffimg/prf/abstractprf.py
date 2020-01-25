# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 20:17:34 2020

@author: fergal
"""

from pdb import set_trace as debug
import numpy as np

#abstractprf.py
class AbstractPrf():
    """All PRF objects must implement a single accessor,
    getPrfForBBox() with the signature below.

    They can implement other accessors, but this is the
    only one that must always be available.
    """

    def __init__(**kwargs):
        pass

    def getPrfForBbox(self, bbox, col0, row0, **kwargs):
        raise NotImplementedError("Abstract base class")


"""
concepts

bbox
    The bounding box of the postage stamp

Model PRF
    An image of what the code expects the PRF to look like for a given
    plate scale

SpecificPrf
    Model image at requested location

ReferencelPrf
    Model PRF precomputed at some reference position on the CCD (e.g at
    the corners)

SubsampledPrf
    An image from which a reference model Prf can be derived

"""

class AbstractLookupPrf():
    """Store and lookup a previously computed PRF function

    This abstract class is created to share as much functionality
    as possible between Kepler and TESS

    To get the recorded prf, use
    getPrfForBbbox(), although this function doesn't work in the base class.
    See docs in that method for more details

    Other functions are in the works to map that PRF onto
    a given mask.

    Todo
    --------
    This class is ripe of optimization with numba
    """

    def __init__(self, path):
        self.path = path
        self.cache = dict()

        self.gridSize = None


    def getPrfForBbox(self, bbox, col, row, *args):
        """Get the prf for an image described by a bounding box

        This function requires as input a function to look up the PRF for a given col row
        (getPrfFunc). This function is not implemented in the base class, as it will be
        specific to each mission. The design intent is that you override getPrfForBbox()
        in the daughter class where
        you define the lookup function, then calls the parent class method.
        See KeplerPrf() for an example

        Input:
        -----------
        col, row
            (floats) Centroid for which to generate prf
        bboxIn
            (int array). Size of image to return. bbox
            is a 4 elt array of [col0, col1, row0, row1]
        getPrfFunc
            (function) Function that returns a PRF object
            This function must have the signature
            ``(np 2d array = getPrfFunc(col, row, *args)``

        Optional Inputs
        -----------------
        Any optional inputs get passed to getPrfFunc


        Returns:
        ----------
        A 2d numpy array of the computed prf at the
        given position.

        Notes:
        ------------
        If bbox is larger than the prf postage stamp returned,
        the missing values will be filled with zeros. If
        the bbox is smaller than the postage stamp, only the
        requestd pixels will be returned
        """

        nColOut = bbox[1] - bbox[0]
        nRowOut = bbox[3] - bbox[2]
        imgOut = np.zeros( (nRowOut, nColOut) )

        #Location of origin of bbox relative to col,row.
        #This is usually zero, but need not be.
        colOffsetOut = (bbox[0] - np.floor(col)).astype(np.int)
        rowOffsetOut = (bbox[2] - np.floor(row)).astype(np.int)

        interpPrf = self.getPrfAtColRow(col, row, *args)
        nRowPrf, nColPrf = interpPrf.shape
        colOffsetPrf = -np.floor(nColPrf/2.).astype(np.int)
        rowOffsetPrf = -np.floor(nRowPrf/2.).astype(np.int)

        di = colOffsetPrf - colOffsetOut
        i0 = max(0, -di)
        i1 = min(nColOut-di , nColPrf)
        if i1 <= i0:
            raise ValueError("Central pixel column not in bounding box")
        i = np.arange(i0, i1)
        assert(np.min(i) >= 0)

        dj = rowOffsetPrf - rowOffsetOut
        j0 = max(0, -dj)
        j1 = min(nRowOut-dj, nRowPrf)
        if j1 <= j0:
            raise ValueError("Central pixel row not in bounding box")
        j = np.arange(j0, j1)
        assert(np.min(j) >= 0)

        #@TODO: figure out how to do this in one step
        for r in j:
            imgOut[r+dj, i+di] = interpPrf[r, i]

        return imgOut


    def interpolatePrf(self, refPrfArray, col, row, evalCols, evalRows):
        """Interpolate between 4 images to find the best PRF at col, row

        This is a private function of the class.

        Inputs
        ------
        refPrfArray
            Array of reference Prfs.
        evalCols
            reference columns for at which the model prfs in refPrfArray
            are computed for.
        evalRows
            reference rows for at which the model prfs in refPrfArray
            are computed for.

        Returns
        -------
        A specific Prf at the location col, row, computed by linear interpolation
        of the reference prfs in two dimensions.

        Note
        --------
        The variable names in this function assume that the refPrfArray
        is ordered as (left-bottom, right-bottom, left-top, right-top),
        but as long as the same order is used for refPrfArray, evalCols,
        and evalRows, the function will work.
        """

        assert len(refPrfArray)== len(evalCols)
        assert len(refPrfArray) == len(evalRows)

        #These asserts are true for Kepler and TESS. May not always
        #be true? They check that reference PRFs are arranged in a square
        #with sides parallel to the sides of the CCD.
        #If these assumptions are not true, the rest of the code does not
        #apply, and needs to be modified
        assert evalCols[0] == evalCols[2]
        assert evalCols[1] == evalCols[3]
        assert evalRows[0] == evalRows[1]
        assert evalRows[2] == evalRows[3]

        p11, p21, p12, p22 = refPrfArray
        c0, c1 = evalCols[:2]
        r0, r1 = evalRows[1:3]

        assert c0 != c1
        assert r0 != r1

        dCol = (col-c0) / (c1-c0)
        dRow = (row-r0) / (r1 - r0)

        #Intpolate across the rows
        tmp1 = p11 + (p21 - p11) * dCol
        tmp2 = p12 + (p22 - p12) * dCol

        #Interpolate across the columns
        out = tmp1 + (tmp2-tmp1) * dRow
        return out
