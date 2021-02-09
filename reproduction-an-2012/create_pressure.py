""" Creates a file for D3D-FM-FLOW with the pressure-field as in the paper An et al. (2012) """

# fix for systems without display
import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp


### Filenames
current_dir = os.path.dirname(os.path.realpath(__file__))


### Parameters
t0 = 10000.  # default: 10000
U = 50.  # default: 50
a = 200000.  # default: 200000
p0 = 2000.  # default: 2000

# cross shore (meters)
x_min = 0.
x_max = 1e6
x_step = 1e4

# along shore (meters)
y_min = -1e7
y_max = 1e7
y_step = x_step

# time (hours)
t_min = 0
t_max = 70 * 3600.
t_step = 1 * 3600.

# command line
parser = argparse.ArgumentParser(
    description="Create pressure field for use with D3D-FM"
)
parser.add_argument(
    "--filename",
    "-f",
    help="Filename (without extension)",
    default="pressure",
    type=str
)
parser.add_argument(
    "--dx",
    help="Grid spacing (in meters)",
    default=None
)
parser.add_argument(
    "--dt",
    help="Time step (in hours)",
    default=None
)

args = parser.parse_args()
filename = str(args.filename)
if args.dx is not None:
    x_step = args.dx
    y_step = x_step
if args.dt is not None:
    t_step = args.dx * 3600.


if filename.endswith(".amp"):
    filename = filename.strip(".amp")
figurename = f"{current_dir}/test_{filename}.jpg"
filename = f"{current_dir}/{filename}.amp"


### Function
def pressure(x, y, t, t0=t0, U=U, a=a, p0=p0):
    return (
        p0
        * (1. - np.exp(-t / t0))
        * np.exp(-(x**2. + (y - U * t)**2. ) / a**2.)
    )


### Compute field
x_num = int((x_max - x_min) / x_step + 1)
y_num = int((y_max - y_min) / y_step + 1)

print("Computing pressure field")
x = np.linspace(x_min, x_max, x_num, dtype=np.float)
y = np.linspace(y_min, y_max, y_num, dtype=np.float)
t = np.arange(t_min, t_max+1, t_step, dtype=np.float)

# loop over all t, y and x
p = np.zeros((t.size, y.size, x.size))
for i in range(t.size):
    for j in range(y.size):
        for k in range(x.size):
            p[i,j,k] = pressure(x[k], y[j], t[i], t0, U, a, p0)

# remove zero-columns and zero-rows
ix = np.where(~ np.all(np.isclose(p, 0), axis=(0,1)))[0]
iy = np.where(~ np.all(np.isclose(p, 0), axis=(0,2)))[0]
x = x[ix]
y = y[iy]
p = p[:,:,ix][:,iy,:]


### Write field
print("Writing pressure field")
data = fp.convert_to_xarray(t, x, y, p)
del t, x, y, p
fp.write_pressure(data, filename)


### Visualise field
print("Plotting pressure field")
fig, ax = plt.subplots(1, 5, sharey=True)
fig.set_tight_layout(True)
im = [None] * 5
for i in range(5):
    idx = data.t.size * i // 5
    im[i] = ax[i].contourf(data.x.values/1000., data.y.values/1000., data.values[idx,:,:], vmin=0, vmax=p0)
    ax[i].set_title(f"$t = {data.t.values[idx] : 0.0f}$s")
    ax[i].set_xlabel("x [km]")
    ax[i].set_xlim([0, 500])
ax[0].set_ylabel("y [km]")
fig.colorbar(im[-2], ax=ax[-1])
plt.savefig(figurename)
plt.show()
