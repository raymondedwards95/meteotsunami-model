
import os
import sys

import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb


###
filename = os.path.dirname(os.path.realpath(__file__)) + "\\bathymetry.xyb"

###
slope = 1./400.

###
def depth(x, alpha=1/400):
    return -1. * alpha * x

###
x = np.linspace(0, 1e6, 5)
y = np.linspace(-1e7, 1e7, 5)
xx, yy = np.meshgrid(x, y)
zz = depth(xx, slope)

###
data = fb.convert_to_xarray(x, y, zz)
fb.write_bathymetry(data, filename)
