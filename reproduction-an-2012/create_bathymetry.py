""" Creates a file for D3D-FM-FLOW with a bathymetry as in paper An et al. (2012) """

import os
import sys

import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb


### Parameters
cases = [0, 36, 41, 42, 43]
slopes = [1./400., 1./800., 0., 0., 0.]
depths = [0, 0, 250, 100, 500]


###
def simple_slope(x, alpha=1/400, d0=0.):
    return -1. * (d0 + alpha * x)


### Paths
bathymetry_dir = os.path.dirname(os.path.realpath(__file__)) + "/bathymetry"
os.makedirs(bathymetry_dir, exist_ok=True)


### Grid
x = np.linspace(0, 1e6, 5)
y = np.linspace(-1e7, 1e7, 5)
xx, yy = np.meshgrid(x, y)


### Computations
for i in range(len(cases)):
    ## Take parameters
    case = cases[i]
    slope = slopes[i]
    depth_0 = depths[i]
    print(f"###\nCreating bathymetry for {case=:02.0f}, with {slope=} and {depth_0=}")

    ## Compute bathymetry
    zz = simple_slope(xx, slope, depth_0)

    ## Write to file
    data = fb.convert_to_xarray(x, y, zz)
    fb.write_bathymetry(data, f"{bathymetry_dir}/repr_{case:02.0f}.xyb")

    ## Visualize
    fb.plot_bathymetry(data, f"{bathymetry_dir}/fig_repr_{case:02.0f}")


# plt.show()
print("Finished creating bathymetry-files")
