""" Additional analysis of waterlevel data """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.optimize import curve_fit

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.utilities as fu


def exp_decay(x, k0, y0):
    return y0 * np.exp(-k0 * x)


def compute_decay_parameter(data, y, t):
    popt, _ = curve_fit(
        exp_decay,  # f
        data["x"],  # xdata
        data["wl"].interp(t=fu.to_timestr(t), y=y),  # ydata
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


def compute_wave_periods(data, y, x=None, crests=True):
    """ Computes the wave period at given x,y """
    t_idx = fu.find_peaks_const_y(data, y, x=x, crests=crests)
    t_idx = np.sort(t_idx)

    if t_idx.size < 2:
        print("No wave period can be computed!")
        return np.nan

    periods = data["t"][t_idx].diff(dim="t").values.astype("timedelta64[s]").astype(float) / 3600

    return periods


def compute_wave_lengths(data, t, x=None, crests=True):
    y_idx = fu.find_peaks_const_t(data, t, x=x, crests=crests)
    y_idx = np.sort(y_idx)

    if y_idx.size < 2:
        print("No wave period can be computed!")
        return np.nan

    lengths = data["y"][y_idx].diff(dim="y").values

    return lengths


def spectral_analysis(data, x, y, t):
    raise NotImplementedError


if __name__ == "__main__":
    _test_compute_decay_parameter(showfig=True)
