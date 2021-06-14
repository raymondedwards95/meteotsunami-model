""" Additional functions """

import datetime

import numpy as np
import xarray as xr


def _filter_peaks(wl, wl_idx, wl_std, window):
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
                if wl[i] > wl_std:  
                    idx.append(i)
    return np.array(idx)


@np.vectorize
def to_timestr(seconds):
    """ Converts time in seconds since reference to a date-string """
    return datetime.datetime.fromtimestamp(seconds).strftime("%Y-%m-%d %H:%M:%S")


def find_peaks_const_y(data, y, x=None):
    """ Finds the times at which the peaks of the waves are present for a given y-coordinate. 
    This can be used to determine the speed at which the waves travel. 

    Input:
        data:   Dataset that contains coordinates and data
        y:      y-coordinate
    
    Parameters:
        x:      x-coordinate

    Output:     indices corresponding to the largest maxima
    """
    window = 10  # number of neighbours to consider

    if x is None:
        x = data["x"].max().values / 50.  # random scaling? close to shore

    if not np.isscalar(y):
        raise ValueError(f"{y=} is not a number")
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")
    if not isinstance(data, (xr.Dataset)):
        raise ValueError(f"{data=} is not a Dataset")
    
    # Extract data
    wl = data["wl"].interp(x=x, y=y)
    wl_std = np.std(wl)

    # Find largest values
    wl_idx = wl.argsort().values[::-1]

    # Find distinct peaks
    return _filter_peaks(wl, wl_idx, wl_std, window)


def find_peaks_const_t(data, t, x=None):
    """ Finds the y-coordinates of the maxima in water level at a given time. 
    This information can be used to compute the wavelength of the wave. 

    Input:
        data:   Dataset that contains coordinates and data
        t:      t-coordinate
    
    Parameters:
        x:      x-coordinate

    Output:     indices corresponding to the largest maxima
    """
    window = 50  # number of neighbours to consider

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
    wl_idx = wl.argsort().values[::-1]

    # Find distinct peaks
    return _filter_peaks(wl, wl_idx, wl_std, window)
