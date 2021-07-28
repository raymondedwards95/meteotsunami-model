""" Additional functions """

import datetime
import os
import sys

import numpy as np
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


def _filter_peaks(wl, wl_idx, wl_std, window, factor):
    """ Find the largest values, remove large values when they are close """
    idx = []  # list of indices corresponding to maxima
    for i in wl_idx:
        if len(idx) == 0:
            # list is empty
            idx.append(i)
        else:
            for j in idx:  
                # check if it is neighbour of previous
                if np.abs(i-j) < window:
                    break
            else:  # i is a new peak
                # check if peak is large enough
                if wl[i] * factor > wl_std:
                    idx.append(i)
    return np.array(idx)


@np.vectorize
def to_timestr(seconds):
    """ Converts time in seconds since reference to a date-string """
    return datetime.datetime.fromtimestamp(seconds).strftime("%Y-%m-%d %H:%M:%S")


def find_peaks_const_y(data, y, x=None, window=10, crests=True, variable="wl"):
    """ Finds the times at which the local maxima (or minima) of the waves are present for a given y-coordinate

    Input:
        data:   Dataset that contains coordinates and data
        y:      y-coordinate
    
    Parameters:
        x:      x-coordinate
        window: expected width of a peak (default: 10)
        factor: find crests (True) or throughs (False)
        variable: variable to use, e.g. "wl", "u" or "v"

    Output:     
        t_idx:  indices corresponding to the time of the largest maxima
    """
    if x is None:
        x = data["x"].max().values / 50.  # random scaling? close to shore

    if not np.isscalar(y):
        raise ValueError(f"{y=} is not a number")
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")
    if not isinstance(data, (xr.Dataset)):
        raise ValueError(f"{data=} is not a Dataset")

    if not variable in ["wl", "u", "v", "p"]:
        raise ValueError(f"Expected {variable=} to be 'wl', 'u', 'v' or 'p'")
    
    # Extract data
    var = data[variable].interp(x=x, y=y)
    var_std = np.std(var)

    # Find largest values
    var_idx = var.argsort().values  # sorts from low to high
    factor = -1

    if crests:
        var_idx = np.flip(var_idx)
        factor = 1

    # Find distinct peaks
    t_idx = _filter_peaks(var, var_idx, var_std, window, factor)
    return t_idx


def find_peaks_const_t(data, t, x=None, window=50, crests=True):
    """ Finds the y-coordinates of the maxima in water level at a given time

    Input:
        data:   Dataset that contains coordinates and data
        t:      t-coordinate
    
    Parameters:
        x:      x-coordinate
        window: expected width of a peak (default: 50)
        crests: find crests (True) or throughs (False)

    Output:     
        y_idx:  indices corresponding to the y-coordinates of the largest maxima
    """
    if x is None:
        x = data["x"].max().values / 50.  # random scaling? close to shore

    if np.isscalar(t):
        t = to_timestr(t)
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")
    if not isinstance(data, (xr.Dataset)):
        raise ValueError(f"{data=} is not a Dataset")
    
    # Extract data
    wl = data["wl"].interp(x=x, t=t)
    wl_std = np.std(wl)

    # Find largest values
    wl_idx = wl.argsort().values  # sorts from low to high
    factor = -1

    if crests:
        wl_idx = np.flip(wl_idx)
        factor = 1

    # Find distinct peaks
    y_idx = _filter_peaks(wl, wl_idx, wl_std, window, factor)
    return y_idx
