""" Creates a bathymetry file for D3D-FM-FLOW """

import os
import sys
import time

import numpy as np

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb
# fmt: on

t0 = time.perf_counter()
print(f"\nStart creating bathymetry-files for exp")


# Function
def exponential_shelf(x, h=20, a=1e-5):
    return -1.0 * h * (1.0 - np.exp(-1.0 * a * x))


# Paths
script_dir = os.path.dirname(os.path.realpath(__file__))
bathymetry_dir = f"{script_dir}/bathymetry"
os.makedirs(bathymetry_dir, exist_ok=True)


# Compute bathymetry
x = np.sort(np.unique(np.concatenate([np.arange(0, 1, 0.5), np.logspace(0, 6, 100)])))
y = np.linspace(-2e6, +3e6, 5)
xx, yy = np.meshgrid(x, y)
zz = exponential_shelf(xx, a=5e-4)


# Write to file
data = fb.convert_to_xarray(x, y, zz, savename=f"{bathymetry_dir}/exp_00.nc")
fb.write_bathymetry(data, f"{bathymetry_dir}/exp_00.xyb")


# Visualize
fb.plot_bathymetry(data, f"{bathymetry_dir}/exp_00", xmax=10e3, keep_open=False)


# End
t1 = time.perf_counter()
print(f"\nFinished creating bathymetry-files for exp in {t1-t0:0.1f} seconds")
