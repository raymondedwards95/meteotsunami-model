""" Creates a file for D3D-FM-FLOW with a bathymetry as in paper An et al. (2012) """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb


###
bathymetry_dir = os.path.dirname(os.path.realpath(__file__)) + "/bathymetry"
os.makedirs(bathymetry_dir, exist_ok=True)
filename = f"{bathymetry_dir}/repr_00.xyb"


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


###
i = y.size // 2
plt.figure()
plt.plot(x, zz[i, :])
plt.title(f"Bottom cross-section at $y={y[i]}$")
plt.savefig(f"{bathymetry_dir}/profile_bath_repr_00.jpg")

plt.figure()
plt.contourf(x, y, zz)
plt.axhline(y[i], color="gray")
plt.title(f"Bottom contours")
plt.savefig(f"{bathymetry_dir}/contour_bath_repr_00.jpg")

# plt.show()
