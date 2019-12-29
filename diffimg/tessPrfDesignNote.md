
## Introduction
The TESS project releases a data file containing realisations of the TESS Pixel Response Function (PRF) at a range of locations across the CCD. The DAVE project uses this PRF file to reconstruct the PRF at arbitrary inter- and intra- pixel locations.

We found that extracting a PRF object from the supplied data file is non-trivial and, given the interest in rapid analysis of the newly released TESS data, we have decided to advertise our Python code for reading and interpolating these files  to the community before rigourous testing is complete.

**This code is made available on a best effort basis**, and we anticipate that some bugs will be found as it gets more use. We strongly suggest you sanity check any results from this code before drawing any conclusions from the results.

If you do use this code, please let us know so we can alert you when updates occur.


## Usage
Download the TESS PRF matlab files from [MAST](https://archive.stsci.edu/tess/all_products.html) . There are 16 files, one for each camera/cdd combination

Download the modules `AbstractPrfLookup.py` and `tessprf.py` from [the DAVE github repo](https://github.com/barentsen/dave/tree/tess-firstlook/diffimg) (we haven't properly packaged DAVE yet)

Create Prf object

    import tessprf
    prfPath  "/path/to/prf/matlab/files"
    ccd, camera, sector = 1, 1, 1
    column, row = 123, 456
    prjObj = tessprf.TessPrf(column, row, ccd, camera, sector)
    
    prfObj.getPrfAtColRow(column, row, ccd, camera, sector)
    
    
## Notes
The PRF is modeled on a grid of 5x5 locations on the CCD, and 9x9 sub-pixel locations at each grid point. To reconstruct the PRF for an arbitrary location on CCD, we should interpolate between grid locations, and intra-pixel locations. Currently the code assumes the difference between the model at adjacence intra-pixel locations is small, and instead returns the PRF model for the nearest intra-pixel locations. If we discover this performance is not sufficient we will remedy this situation


## Performance
We tested the behaviour of our class for 3 locations at the bottom left corner of ccd-camera = 1-1. Astigmatism is a problem in this part of the field of view, especially low values of column and row. This changing astigmatism gives us an opportunity to test the peformance of our class.

We show the results in the following figures. In each figure the left panel shows the data from a single FFI (tess2018218052942-s0001-1-1-0120-s_ffic.fits, available from MAST) in units of data number, the middle panel shows our model at that location, and the rightmost panel shows the difference image. The score is defined as the average absolute difference between the image and the model in data number

![ ](./tess1-1-1-eg1.png  "Figure 1")
Figure 1: Data, model, and residuals for star at (col, row) on CCD 1 of Camera 1
The data and model are shown in a log scaling, the residuals in a linear one. The model is not perfect, but the residuals are a small fraction of the original flux. The shape of the PRF is due to the large astigmatism this far off axis on the focal plane.

![ ](./tests/tess1-1-1-eg2.png  "Figure 2")
Figure 2: Same as Figure 1, but for a PRF at (col, row).

![ ](./tests/tess1-1-1-eg3.png  "Figure 3")
Figure 3: Same as Figure 1, but for a star at (col, row), much closer to the centre of the camera.