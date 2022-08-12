""" Creates a bathymetry file for D3D-FM-FLOW """

import os
import sys

import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb


### Function
def exponential_shelf(x, h=20, a=1e-5):
    return -1. * h * (1. - np.exp(- 1. * a * x))


### Paths
bathymetry_dir = os.path.dirname(os.path.realpath(__file__)) + "/bathymetry"
os.makedirs(bathymetry_dir, exist_ok=True)


### Compute bathymetry
x = np.concatenate(([0], np.logspace(0, 6, 101)))
y = np.linspace(-2e6, +3e6, 5)
xx, yy = np.meshgrid(x, y)
zz = exponential_shelf(xx, a=5e-4)


### Write to file
data = fb.convert_to_xarray(x, y, zz)
fb.write_bathymetry(data, f"{bathymetry_dir}/exp_00.xyb")


### Visualize
fb.plot_bathymetry(data, f"{bathymetry_dir}/fig_exp_00", xmax=10e3)

print("Finished creating bathymetry-files")
