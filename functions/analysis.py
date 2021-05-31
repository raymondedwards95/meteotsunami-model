""" Additional analysis of waterlevel data """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.optimize import curve_fit

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.visualisation as fv


def exp_decay(x, k0, y0):
    return y0 * np.exp(-k0 * x)


def compute_decay_parameter(data, y, t):
    popt, _ = curve_fit(
        exp_decay,  # f
        data["x"],  # xdata
        data["wl"].interp(t=fv.to_timestr(t), y=y),  # ydata
        p0=[1./100000., 1.],  # p0
        bounds=(0., [10., 10.])  # bounds (min, max)
    )
    k0, y0 = popt
    return k0, y0


def _test_compute_decay_parameter(showfig=True):
    k0 = 1./12.
    y0 = 0.7

    x = np.arange(100)
    y = exp_decay(x, k0, y0) + np.random.normal(scale=0.01, size=x.size)

    (k0_est, y0_est), _ = curve_fit(
        exp_decay,
        x,
        y
    )

    plt.figure()
    plt.scatter(x, y)
    plt.plot(x, exp_decay(x, k0_est, y0_est), "C1", label=f"$1/k_0={1/k0_est:0.1f}$; $y_0={y0_est:0.2f}$")
    plt.legend()
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.show()


if __name__ == "__main__":
    _test_compute_decay_parameter(showfig=True)
