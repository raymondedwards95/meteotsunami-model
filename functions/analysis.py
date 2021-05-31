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

