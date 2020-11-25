
import os
import sys

import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp


###
filename = os.path.dirname(os.path.realpath(__file__)) + "\\pressure.amp"

###
t0 = 10000.
U = 50.
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
x = np.linspace(0, 1e6, 51)
y = np.linspace(-1e7, 1e7, 75)
t = np.arange(0, 70, 1) * 3600.

p = np.zeros((t.size, y.size, x.size))
for i in range(t.size):
    for j in range(y.size):
        for k in range(x.size):
            p[i,j,k] = pressure(x[k], y[j], t[i])

###
data = fp.convert_to_xarray(t, x, y, p)
fp.write_pressure(data, filename)
