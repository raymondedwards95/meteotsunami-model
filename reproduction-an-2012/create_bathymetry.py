""" Creates a file for D3D-FM-FLOW with a bathymetry as in paper An et al. (2012) """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb


###
cases = [0, 36]
slopes = [1./400., 1./800.]


###
def depth(x, alpha=1/400):
    return -1. * alpha * x


###
bathymetry_dir = os.path.dirname(os.path.realpath(__file__)) + "/bathymetry"
os.makedirs(bathymetry_dir, exist_ok=True)


###
x = np.linspace(0, 1e6, 5)
y = np.linspace(-1e7, 1e7, 5)
xx, yy = np.meshgrid(x, y)


###
for i in range(len(cases)):
    case = cases[i]
    slope = slopes[i]
    print(f"###\nCreating bathymetry for case {case:02.0f}")

    zz = depth(xx, slope)
    data = fb.convert_to_xarray(x, y, zz)
    fb.write_bathymetry(data, f"{bathymetry_dir}/repr_{case:02.0f}.xyb")
    fb.plot_bathymetry(data, f"{bathymetry_dir}/fig_repr_{case:02.0f}")


# plt.show()
print("Finished creating bathymetry-files")
