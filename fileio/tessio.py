
"""
Created on Tue Nov 27 20:37:41 2018

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from dave.fileio.AbstractMast import MastArchive
import os

class TessDvtLocalArchive(MastArchive):
    """
    Shim class to load DV Timeseries data from a local disk
    
    /Users/smullally/TESS/TCEs/Sector2/dvt/tess2018235142541-s0002-s0002-0000000100100827-00109_dvt.fits”

    tess2018206190142-s0001-s0001-0000000471013508-00106_dvt.fits        

    """    
    def __init__(self, path):
        MastArchive.__init__(self, "/dev/null", 'www.example.com')

        self.path = path

        self.filePrefix = dict()
        self.filePrefix[1] = "tess2018206190142-s0001-s0001"        
        self.filePrefix[2] = "tess2018235142541-s0002-s0002"
        
        self.fileSuffix = dict()
        self.fileSuffix[1] = "00106_dvt.fits"
        self.fileSuffix[2] = "00109_dvt.fits"
        
    def getDvt(self, sector, tic, *args, **kwargs):
        
        prefix= self.filePrefix[sector]
        suffix = self.fileSuffix[sector]
        localUrl = "%s-%016i-%s" %(prefix, int(tic), suffix)
        localUrl = os.path.join(self.path, localUrl)

        if not os.path.exists(localUrl):
            raise IOError("File not found: %s" %(localUrl))
            
        return self.parse(localUrl, *args, **kwargs)
            
