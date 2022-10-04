""" Helper functions for analysis and visualisation

Main functions:
    to_timestr
    find_peaks_const_y
    find_peaks_const_t
"""

import datetime
import os
import sys

import numpy as np
import numpy.typing as npt
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
# fmt: on


def _filter_peaks(
    wl: npt.ArrayLike,
    wl_idx: npt.ArrayLike,
    wl_std: Numeric,
    window: Integer,
    factor: Numeric,
) -> np.ndarray:
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
def to_timestr(seconds: Numeric) -> str:
    """ Converts time in seconds since reference to a date-string """
    return datetime.datetime.fromtimestamp(seconds).strftime("%Y-%m-%d %H:%M:%S")


def find_peaks_const_y(
    data: xr.Dataset,
    y: Numeric,
    x: Numeric = None,
    window: Integer = 10,
    crests: bool = True,
    variable: str = "wl",
) -> np.ndarray:
    """ Finds the times at which the local maxima (or minima) of the waves are present for a given y-coordinate

    Input:
        `data`:     Dataset that contains coordinates and data
        `y`:        y-coordinate
        `x`:        x-coordinate

    Options:
        `window`:   expected width of a peak (default: 10)
        `factor`:   find crests (True) or throughs (False)
        `variable`: variable to use, e.g. "wl", "u" or "v"

    Output:
        `t_idx`:    indices corresponding to the time of the largest maxima
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

    # End
    return t_idx


def find_peaks_const_t(
    data: xr.Dataset,
    t: Numeric,
    x: Numeric = None,
    window: Integer = 50,
    crests: bool = True,
) -> np.ndarray:
    """ Finds the y-coordinates of the maxima in water level at a given time

    Input:
        `data`:     Dataset that contains coordinates and data
        `t`:        t-coordinate
        `x`:        x-coordinate

    Options:
        `window`:   expected width of a peak (default: 50)
        `crests`:   find crests (True) or throughs (False)

    Output:
        `y_idx`:    indices corresponding to the y-coordinates of the largest maxima
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


def find_local_maxima_y(
    dataset: xr.Dataset,
    t: Numeric,
    x: Numeric,
    variable: str = "wl",
    minima: bool = False,
) -> np.ndarray:
    """ Finds y-coordinates of the local maxima for fixed (t,x)

    Input:
        `dataset`:  dataset containing gridded model output
        `t`:        t-coordinate
        `x`:        x-coordinate

    Options:
        `variable`: name of variable
        `minima`:   find local minima instead of maxima
    """
    # Extract data
    data = dataset[variable].interp(x=x, t=t)

    # Get lower limit
    data_std = np.std(data)

    # Invert data if minima are asked
    if minima:
        data = -1. * data

    # Find slope of data
    data_diff = data \
        .differentiate("y") \
        .data  # data converts from xarray to a plain numpy array

    # Find indices of zero-crossings of derivative, i.e. local minima/maxima of data, by making use of sign-changes
    y_idx = np.where(
        (data_diff[:-1] * data_diff[1:]) < 0
    )[0]

    # Find indices for waves larger than a certain lower limit
    return y_idx[data[y_idx] > data_std]


if __name__ == "__main__":
    # Define paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    anim_dir = f"{script_dir}/tests/anim"
    os.makedirs(anim_dir, exist_ok=True)

    # Read data
    data = xr.open_dataset(
        f"{script_dir}/../reproduction-an-2012/output/data_repr_00.nc")

    # Test functions
    find_local_maxima_y(
        dataset=data,
        t=to_timestr(3600.*42.),
        x=1e4,
    )
