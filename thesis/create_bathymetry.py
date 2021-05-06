""" Creates a bathymetry file for D3D-FM-FLOW """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb


###
@np.vectorize
def exponential_shelf(x, h=20, a=1e-5):
    return -1. * h * (1. - np.exp(- 1. * a * x))


###
bathymetry_dir = os.path.dirname(os.path.realpath(__file__)) + "/bathymetry"
os.makedirs(bathymetry_dir, exist_ok=True)
filename = f"{bathymetry_dir}/exp_00.xyb"


###
x = np.linspace(0, 1e6, 51)
y = np.linspace(-1e7, 1e7, 5)
xx, yy = np.meshgrid(x, y)
zz = exponential_shelf(xx)


###
data = fb.convert_to_xarray(x, y, zz)
fb.write_bathymetry(data, filename)


###
i = y.size // 2
plt.figure()
plt.plot(x, zz[i, :])
plt.title(f"Bottom cross-section at $y={y[i]}$")
plt.savefig(f"{bathymetry_dir}/profile_bath_exp_00.jpg")

plt.figure()
plt.contourf(x, y, zz)
plt.axhline(y[i], color="gray")
plt.title(f"Bottom contours")
plt.savefig(f"{bathymetry_dir}/contour_bath_exp_00.jpg")

# plt.show()
