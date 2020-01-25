
import numpy as np
import astropy.io.fits as pyfits
import os

__version__ = "$Id: prf.py 1964 2015-02-27 18:33:11Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/prf.py $"


from AbstractLookupPrf import AbstractLookupPrf


class KeplerPrf(AbstractLookupPrf):
    """Return the expected PRF for a point source in the Kepler field
    based on mod out, and centroid.

    Note:
    --------
    For speed, this function caches some data for each module to
    save reading it from file for each call. If you request the PRF
    from many mod-outs you may find your memory increases dramatically.

    """

    def __init__(self, path, mod, out):
        AbstractLookupPrf.__init__(self, path)
        self.gridSize = 50

        assert (1<= mod) and (mod <= 24)
        assert (1<= out) and (out <= 4)
        self.mod = mod
        self.out = out

    def getPrfForBbox(self, bbox, col, row, mod, out):
        """Get PRF for a bounding box.

        See documentation in the same method in the parent class
        """
        self.validateInputs(col, row, mod, out)
        args = [mod, out]

        return AbstractLookupPrf.abstractGetPrfForBbox(bbox, col, row, *args)

    def validateInputs(self, col, row):
        return True
        #@TODO Check col, row in range


    def getPrfAtColRow(self, col, row):
        """Compute the model prf for a given module, output, column row

        This is the workhorse function of the class. For a given mod/out,
        loads the subsampled PRFs at each corner of the mod-out, extracts
        the appropriate image for the given subpixel position, then interpolates
        those 4 images to the requested column and row.

        Note:
        --------
        For speed, this function caches some data for each module to
        save reading it from file for each call. If you request the PRF
        from many mod-outs you may find your memory increases dramatically.

        """

        mod, out = self.mod, self.out

        #Load subsampled PRFs from cache, or from file if not previously read in
        key = "%02i-%02i" %(mod, out)
        if key not in self.cache:
            self.cache[key] = self.getSubSampledPrfs(mod, out)

        fullPrfArray = self.cache[key]
        regPrfArray = self.getRegularlySampledPrfs(fullPrfArray, col,row)
        bestRegPrf = self.interpolateRegularlySampledPrf(regPrfArray, \
            col, row)

        return bestRegPrf


    def getRegularlySampledPrfs(self, fullPrfArray, col, row):
        regArr = []
        for i in range(5):
            tmp = self.getSingleRegularlySampledPrf(fullPrfArray[i], \
                col, row)
            regArr.append(tmp)
        return regArr


    def getSingleRegularlySampledPrf(self, singleFullPrf, col, row):
        """Extract out an 11x11* PRF for a given column and row for
        a single 550x550 representation of the PRF

        Note documentation for this function is Kepler specific

        Doesn't interpolate across prfs just takes
        care of the intrapixel variation stuff

        Steve stored the prf in a funny way. 50 different samplings
        of an 11x11 pixel prf*, each separated by .02 pixels are
        stored in a 550x550 grid. Elements in the grid 50 spaces
        apart each refer to the same prf.


        *Sometimes a 15x15 pixel grid

        Notes:
        ---------
        None of the details of how to write this function are availabe
        in the external documents. It was all figured out by trial and error.
        """
        gridSize = self.gridSize

        #The sub-pixel image from (0.00, 0.00) is stored at x[49,49].
        #Go figure.
        colIndex = ((1-np.remainder(col, 1)) * gridSize).astype(np.int32)-1
        rowIndex = ((1-np.remainder(row, 1)) * gridSize).astype(np.int32)-1

        nColOut, nRowOut = singleFullPrf.shape
        nColOut /= float(gridSize)
        nRowOut /= float(gridSize)

        iCol = colIndex + (np.arange(nColOut)*gridSize).astype(np.int)
        iRow = rowIndex + (np.arange(nRowOut)*gridSize).astype(np.int)

        #Don't understand why this must be a twoliner
        tmp = singleFullPrf[iRow, :]
        return tmp[:,iCol]


    def interpolateRegularlySampledPrf(self, regPrfArray, col, row):

        #See page 2 of PRF_Description_MAST.pdf for the storage
        #order
        p11, p12, p21, p22, pMid = regPrfArray

        #See instrument hand bok, page 49, Section 4.5
        nCol, nRow = 1099, 1023
        col = int(col) / float(nCol)
        row = int(row) / float(nRow)

        #Intpolate across the rows
        tmp1 = p11 + (p21 - p11) * col
        tmp2 = p12 + (p22 - p12) * col

        #Interpolate across the columns
        out = tmp1 + (tmp2-tmp1)*row
        return out


    def mapPrfToImg(self, bestRegPrf, imgSizeRC, imgOffsetCR, \
            centroidCR):
        """Place a rectangular apeture over the prf

        Note this will fail if the image aperture doesn't wholly
        overlap with bestRegPrf
        """

        #Location of brightest pixel in PRF img. Usually 5,5,
        #but sometimes 7,7
        midPrfRow, midPrfCol = np.floor(np.array(bestRegPrf.shape) / 2.)
        nRowImg, nColImg = imgSizeRC

        #Map midPrf so it lies on the same pixel within the image
        #centroidCR does
        xy = np.array(centroidCR) - np.array(imgOffsetCR)

        c1 = midPrfCol-xy[0] + 1
        c1 = np.arange(c1, c1+imgSizeRC[1], dtype=np.int32)

        r1 = midPrfRow-xy[1] + 1
        r1 = np.arange(r1, r1+imgSizeRC[0], dtype=np.int32)

        tmp = bestRegPrf[r1, :]
        return tmp[:,c1]


    def getSubSampledPrfs(self, mod, out):
        fullPrfArray = []
        for i in range(1,6):
            tmp = self.readPrfFile(mod, out, i)
            fullPrfArray.append(tmp)
        return fullPrfArray


    def readPrfFile(self, mod, out, position):
        """

        position    (int) [1..5]
        """

        filename = "kplr%02i.%1i_2011265_prf.fits" %(mod, out)
        filename = os.path.join(self.path, filename)

        img = pyfits.getdata(filename, ext=position)
        return img


