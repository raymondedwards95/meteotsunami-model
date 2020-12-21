""""""

# fix for systems without display
import matplotlib
matplotlib.use("Agg")

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp


### Filenames
current_dir = os.path.dirname(os.path.realpath(__file__))
filename = current_dir + "/pressure.amp"

### Parameters
t0 = 10000.  # default: 10000
U = 50.  # default: 50
a = 200000.  # default: 200000
p0 = 2000.  # default: 2000

### Function
def pressure(x, y, t, t0=t0, U=U, a=a, p0=p0):
    return (
        p0
        * (1. - np.exp(-t / t0))
        * np.exp(-(x**2. + (y - U * t)**2. ) / a**2.)
    )

### Compute field
print("Computing pressure field")
x = np.linspace(0, 1e6, 101, dtype=np.float)
y = np.linspace(-1e7, 1e7, 501, dtype=np.float)
t = np.arange(0, 70, 1, dtype=np.float) * 3600.

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
plt.savefig(current_dir + "/test_pressure.jpg")
plt.show()
