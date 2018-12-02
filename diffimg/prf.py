
import numpy as np
import astropy.io.fits as pyfits
import os

__version__ = "$Id: prf.py 1964 2015-02-27 18:33:11Z fmullall $"
__URL__ = "$URL: svn+ssh://flux/home/fmullall/svn/kepler/py/prf.py $"

npmap = lambda f, x: np.array(map(f, x))

class AbstractBasePrfLookup(object):
    """Store and lookup a previously computed PRF function

    This abstract class is created in the hope that much of the functionality can
    be reused for TESS.

    To get the recorded prf, use
    getPrfForBbbox(), although this function doesn't work in the base class. See docs
    in that method for more details

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


    def abstractGetPrfForBbox(self, col, row, bboxIn, getPrfFunc, *args):
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

        bbox = np.array(bboxIn).astype(int)
        nColOut = bbox[1] - bbox[0]
        nRowOut = bbox[3] - bbox[2]
        imgOut = np.zeros( (nRowOut, nColOut) )

        #Location of origin of bbox relative to col,row.
        #This is usually zero, but need not be.
        colOffsetOut = (bbox[0] - np.floor(col)).astype(np.int)
        rowOffsetOut = (bbox[2] - np.floor(row)).astype(np.int)

        interpPrf = getPrfFunc(col, row, *args)
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






class KeplerPrf(AbstractBasePrfLookup):
    """Return the expected PRF for a point source in the Kepler field
    based on mod out, and centroid.

    Note:
    --------
    For speed, this function caches some data for each module to
    save reading it from file for each call. If you request the PRF
    from many mod-outs you may find your memory increases dramatically.

    """

    def __init__(self, path):
        AbstractBasePrfLookup.__init__(self, path)
        self.gridSize = 50

    def getPrfForBbox(self, mod, out, col, row, bboxIn):
        """Get PRF for a bounding box.

        See documentation in the same method in the parent class
        """
        args = [mod, out]

        return self.abstractGetPrfForBbox(col, row, bboxIn, self.getPrfAtColRow, *args)


    def getPrfAtColRow(self, col, row, mod, out):
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


from pdb import set_trace as debug

import scipy.io as spio
class TessPrf(AbstractBasePrfLookup):
    
    """Interpolate a TESS PRF image
    
    Caution: This code is under development and is expected to be buggy
    
    TODO: Interplolate, then pull out the regular sampled PRF
    
    Cache interplolation results for speed
    """
    
    def __init__(self, path):
        AbstractBasePrfLookup.__init__(self, path)
        self.gridSize = 9

    def getPrfForBbox(self, sector, camera, ccd, col, row, bboxIn):
        """Get PRF for a bounding box.

        See documentation in the same method in the parent class
        """
        args = [ccd, camera, sector]

        return self.abstractGetPrfForBbox(col, row, bboxIn, self.getPrfAtColRow, *args)

    
    def checkOutOfBounds(self, col, row):
        if col < 45 or col > 2091:
            raise ValueError("Requested column (%i) not on phyiscal CCD [45,2091]" %(col))
            
        if row < 1 or row > 2047:
            raise ValueError("Requested row (%i) not on phyiscal CCD [0,2047]" %(row))

    def getPrfAtColRow(self, col, row, ccd, camera, sector):

        self.checkOutOfBounds(col, row)        
        #Currently, the same PRF model applies to all sectors, so
        #we pin the sector number. If the PRF is recalculated at a later
        #date we'll need some logic here.
        sector = 1
        key = "%1i-%1i-%02i" %(ccd, camera, sector)
        
        if key not in self.cache:
            self.cache[key] = self.readPrfFile(ccd, camera, sector)

        prfObj = self.cache[key]
        prfArray, evalCols, evalRows = self.getRegularlySampledBracketingPrfs(prfObj, col,row)
        bestPrf = self.interpolatePrf(prfArray, \
            col, row, evalCols, evalRows)

#        regPrf = self.getSingleRegularlySampledPrf(bestPrf, col, row)
        return bestPrf

    
    def getRegularlySampledBracketingPrfs(self, prfObj, col, row):
        #Get columns and rows at which PRF was evaluated at
        cr = np.array((col, row))

        nEval = np.sqrt(len(prfObj))
        assert nEval == int(nEval), "PRF grid is not square"
        nEval = int(nEval)
        
        #Read cols and rows at which PRF is evaluated
        evalCol = npmap(lambda x: x.ccdColumn, prfObj)
        evalRow = npmap(lambda x: x.ccdRow, prfObj)
        
        #Sort prfArray so its the right shape for getBracketingIndices()
        evalColrow = np.vstack((evalCol, evalRow)).transpose()
        srt = np.lexsort((evalRow, evalCol))
        evalColRow = evalColrow[srt].reshape((nEval, nEval, 2))
        prfArray = prfObj[srt].reshape((nEval, nEval))
        
        whBracket = getBracketingIndices(evalColRow, cr)        
        
        c0, r0 = [], []
        regPrfArr = []
        for itr in range(4):
            i, j = whBracket[itr]
            #Store column and row
            c0.append( evalColRow[ i, j, 0 ] )
            r0.append( evalColRow[ i, j, 1 ] )

            #Check I did all the book keeping correctly
            assert c0[itr] == prfArray[i,j].ccdColumn    
            assert r0[itr] == prfArray[i,j].ccdRow    
            
            #Pull out the 13x13 image for this location
            regPrf = self.getSingleRegularlySampledPrf(prfArray[i,j], col, row)
            regPrfArr.append(regPrf)

        #More checks: check the order of the locations is correct
        assert c0[0] == c0[2]
        assert c0[1] == c0[3]
        assert r0[0] == r0[1]
        assert r0[2] == r0[3]
        return np.array(regPrfArr), np.array(c0), np.array(r0)



    def getSingleRegularlySampledPrf(self, singlePrfObj, col, row):
        """
        
        Todo
        --------
        This function returns the PRF at the closest point of evaluation. It 
        really should interpolate between points to get a PRF that varies
        more smoothly with intrapixel location.
        """
        #Get the sub-pixel positions at which the PRF is evaluated
        #The sub pixel positions repeat every 9 terms, so we only take the first repitition
        #np.remainder strips out the integer portion.
        #[:gridSize] trims to just the first repitition of the sequence
        gridSize = self.gridSize
        evalCol = np.remainder(singlePrfObj.prfColumn[:gridSize], 1)
        evalRow = np.remainder(singlePrfObj.prfRow[:gridSize], 1)

        #col, rowFrac represent the intraPixel location of the centroid.
        #For example, if (col, row) = (123,4, 987.6), (colFrac, rowFrac) = (.4, .6)
        #PRF evaluated from -.5 to +.5 in the pixel, so I think these half integer
        #offsets are necessary
        colFrac = np.remainder(col, 1) - .5
        rowFrac = np.remainder(row, 1) - .5
        
        colOffset = np.argmin( evalCol - colFrac)
        rowOffset = np.argmin( evalRow - rowFrac)

        #TODO Reference documentation about samplesPerPixel        
        fullImage = singlePrfObj.values/ float(singlePrfObj.samplesPerPixel)
        
        #Number of pixels in regularly sampled PRF. Typically 13x13
        nColOut, nRowOut = fullImage.shape
        nColOut /= float(gridSize)
        nRowOut /= float(gridSize)

        iCol = colOffset + (np.arange(nColOut)*gridSize).astype(np.int)
        iRow = rowOffset + (np.arange(nRowOut)*gridSize).astype(np.int)

        #Don't understand why this must be a twoliner
        tmp = fullImage[iRow, :]
        return tmp[:,iCol]
    
    
    def interpolatePrf(self, regPrfArray, col, row, evalCols, evalRows):
        #Sort inputs so evalCol and evalRows are both monotonically increasing
        #np.interplated needs this to be true

        #TODO Find a faster way to do this
        index = np.array(map( lambda x, y: "%04i-%04i" %(x,y), evalCols, evalRows))   
        
        srt = np.argsort(index) 
        regPrfArray = regPrfArray[srt]
        evalCols = evalCols[srt]
        evalRows = evalRows[srt]

        p11, p12, p21, p22 = regPrfArray
        c0, c1 = evalCols[1:3]
        r0, r1 = evalRows[:2]
        
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
        
    

            
    def readPrfFile(self, ccd, camera, sector):
        
        if camera != 1:
            raise ValueError("Only camera 1 currently available")
            
        fn = "tess2018243163600-00072_035-%i-%i-characterized-prf.mat" %(ccd, camera)
        path = os.path.join(self.path, fn)

        obj = spio.matlab.loadmat(path, struct_as_record=False, squeeze_me=True)
        prfObj = obj['prfStruct']
        return prfObj



def getBracketingIndices(evalColRow, cr):
    """
    Get the indices of `evalColRow` that bracket `cr`

    This is a special function used by TessPrf
    
    This function encapsulates some fairly knotty bookkeeping. Unless something
    is broken you probably want to leave this function well alone
    
    Inputs
    --------
    evalColRow
        (3d np array) See discussion below
    cr
        (2 element np array) The column and row to be bracketed

    Returns
    ----------
    A 4x2 numpy array. Each row represents the indices into
    `evalColRow[,,:]` representing the 4 points in `evalColRow`
    that bracket the location represented by evalColRow
    
    
    Note
    -----
    The model prf is evaluated on a regular grid across the CCD. Each
    grid point can be represented in two coordinate systems; the 
    CCD pixel coordinates (this PRF is evaluated at col,row=24,36,
    and a grid Index (this is the second grid location in column, and
    the third in row). `evalColRow` encodes a mapping from one coord sys
    to the other.
    
    The zeroth dimension of `evalColRow` encodes the column of the grid 
    location (e.g. 2 in the example above). The first dimension
    encodes row of the grid location (3 in the example), the second
    dimension encodes whether the value represents CCD column 
    (`evalColRow[:,:,0]`) or CCD row (`evalColRow[:,:,1]`). The
    value at each array element represents the CCD position (either
    column or row).
    
    The return value of this function is a list of the 4 grid locations
    that bracket the input `cr` in column and row (below left, below right,
    above left, above right)
    
    Example
    ---------
    `evalColRow` consists of 4 points at which the model prf is evalulated
    
    .. code-block:: python
    
        a[0,0,0] =  45
        a[0,0,1] =   1   #Zeroth prf evalulated at (col, row) = (45,1)
        a[0,1,0] =  45
        a[0,1,1] = 128

        a[1,0,0] = 183
        a[1,0,1] =   1
        a[1,1,0] = 183
        a[1,1,1] = 128

        cr = (45, 45)  #Somewhere in the middle
    
    The return value is
    
    .. code-block:: python
    
        [ [0,0], [1,0], [1,0], [1,1] ]
        
    Because these are the indices that bracket the input col,row
    """
    tmp = (evalColRow - cr)
    dist = np.hypot(tmp[:,:,0], tmp[:,:,1])
    wh = np.unravel_index( np.argmin(dist), dist.shape)

    nearestEval = evalColRow[wh]
    delta = cr - nearestEval 

    #Find the 3 other evaluations of the PRF that bracket (col, row)
    tmp = []
    if delta[0] >= 0 and delta[1] >= 0:        #A
        tmp.append( wh )
        tmp.append( wh + np.array((+1, +0)) )
        tmp.append( wh + np.array((+0, +1)) )
        tmp.append( wh + np.array((+1, +1)) )
        
    elif delta[0] < 0 and delta[1] >= 0:       #S
        tmp.append( wh + np.array((-1, +0)) )
        tmp.append( wh )
        tmp.append( wh + np.array((-1, +1)) )
        tmp.append( wh + np.array((+0, +1)) )

    elif delta[0] < 0 and delta[1] < 0:        #T
        tmp.append( wh + np.array((-1, -1)) )
        tmp.append( wh + np.array((-0, -1)) )
        tmp.append( wh + np.array((-1, +0)) )
        tmp.append( wh )

    else:                                      #C
        tmp.append( wh + np.array((-0, -1)) )
        tmp.append( wh + np.array((+1, -1)) )
        tmp.append( wh )
        tmp.append( wh + np.array((+1, +0)) )
    tmp = np.array(tmp)

    #Check the order of values is correct
    c0 = tmp[:,0]
    r0 = tmp[:,1]
    assert c0[0] == c0[2]
    assert c0[1] == c0[3]
    assert r0[0] == r0[1]
    assert r0[2] == r0[3]
    
    #Bounds checking
    assert np.min(tmp) >= 0
    assert np.max(tmp[:,0]) < evalColRow.shape[0]
    assert np.max(tmp[:,1]) < evalColRow.shape[1]

    return tmp

