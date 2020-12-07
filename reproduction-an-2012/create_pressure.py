
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp


###
filename = os.path.dirname(os.path.realpath(__file__)) + "\\pressure.amp"

###
t0 = 10000.
U = - 50.
a = 200000.
p0 = 2000.

###
def pressure(x, y, t, t0=10000, U=50, a=200000, p0=2000):
    return (
        p0
        * (1. - np.exp(-t / t0))
        * np.exp(-(x**2. + (y - U * t)**2. ) / a**2.)
    )

###
print("Computing pressure field")
x = np.linspace(0, 1e6, 51)
y = np.linspace(-1e7, 1e7, 501)
t = np.arange(0, 70, 1) * 3600.

p = np.zeros((t.size, y.size, x.size))
for i in range(t.size):
    for j in range(y.size):
        for k in range(x.size):
            p[i,j,k] = pressure(x[k], y[j], t[i])

###
print("Writing pressure field")
data = fp.convert_to_xarray(t, x, y, p)
fp.write_pressure(data, filename)

###
print("Plotting pressure field")
fig, ax = plt.subplots(1, 5, sharey=True)
fig.set_tight_layout(True)
im = [None] * 5
for i in range(5):
    idx = t.size * i // 5
    im[i] = ax[i].contourf(x/1000., y/1000., p[idx,:,:], vmin=0, vmax=p0)
    ax[i].set_title(f"$t = {t[idx] : 0.0f}$s")
    ax[i].set_xlabel("x [km]")
    ax[i].set_xlim([0, 500])
ax[0].set_ylabel("y [km]")
fig.colorbar(im[-2], ax=ax[-1])
plt.savefig("./example_pressure.jpg")
plt.show()
