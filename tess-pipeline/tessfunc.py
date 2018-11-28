"""
Created on Tue Nov 27 21:44:04 2018

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


from dave.fileio.tessio import TessDvtLocalArchive


def serve(sector, tic, planetNum, localPath):
    ar = TessDvtLocalArchive(localPath)
    
    dvt, hdr = ar.getDvt(tic, sector, ext=planetNum, header=True)

    return dvt, hdr
    