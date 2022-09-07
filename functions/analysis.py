""" Additional analysis of waterlevel and water velocity data

Main functions:
    compute_decay_parameter
    compute_wave_periods
    compute_wave_lengths
    spectral_analysis_1d
    spectral_analysis_2d
"""

import os
import sys
from typing import List, Tuple, Union

import numpy as np
import numpy.typing as npt
import scipy.constants
import xarray as xr
from scipy.optimize import curve_fit

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.utilities as fu
# fmt: on


def dispersion_relation(k: Union[npt.ArrayLike, Numeric], n: Integer = 0, alpha: Numeric = 1/400) -> np.ndarray | Numeric:
    """ Returns the frequency related to the wavenumber
    for shallow water waves on a beach of linear slope [Eckart, 1951]

    Input:
        k:      wavenumbers

    Options:
        n:      mode
        alpha:  slope of beach

    Output:
        f:      frequencies corresponding to wavenumbers
    """
    g = scipy.constants.g
    assert np.isclose(g, 9.81, rtol=1e-2)
    return np.sqrt((g * np.abs(k) / 2. / np.pi) * (2. * n + 1.) * np.tan(alpha))


def exp_decay(x: npt.ArrayLike, k0: Numeric, y0: Numeric) -> np.ndarray:
    """ Returns an exponential decay in x with decay parameter k0 and scale y0

    Input:
        x           horizontal coordinates
        k0          decay parameter
        y0          vertical scale

    Output:
        profile:    profile of the exponential decay at coordinates x
    """
    return y0 * np.exp(-k0 * x)


def compute_decay_parameter(data: xr.Dataset, y: Numeric, t: Numeric) -> Tuple[Numeric]:
    """ Computes the decay parameter and scaling parameter for data at given y and t

    Input:
        data:   dataset containing water level data
        y:      y-coordinate
        t:      t-coordinate in seconds
    """
    wl = data["wl"].interp(t=fu.to_timestr(t), y=y)
    if np.all(np.isnan(wl)):
        return np.nan, 0

    popt, _ = curve_fit(
        exp_decay,  # f
        data["x"],  # xdata
        wl,  # ydata
        p0=[1./100000., 1.],  # p0
    )
    k0, y0 = popt
    return (k0, y0)


def compute_wave_periods(data: xr.Dataset, y: Numeric, x: Numeric = None, crests: bool = True, no_result: Numeric = np.nan) -> List[Numeric]:
    """ Computes the wave period at given x,y

    Input:
        data:       dataset containing all data
        y:          y-coordinate
        x:          x-coordinate

    Options:
        crests:     find crests (True) or throughs (default: False)
        no_result:  value to return if there is no result (default: np.nan)

    Output:
        periods:    list of time-intervals between wave crests
    """
    t_idx = fu.find_peaks_const_y(data, y, x=x, crests=crests)
    t_idx = np.sort(t_idx)

    if t_idx.size < 2:
        print("No wave period can be computed!")
        return np.array([no_result])

    periods = \
        data["t"][t_idx] \
        .diff(dim="t") \
        .values \
        .astype("timedelta64[s]") \
        .astype(float) / 3600.

    return periods


def compute_wave_lengths(data: xr.Dataset, t: Numeric, x: Numeric = None, crests: bool = True, no_result: Numeric = np.nan) -> List[Numeric]:
    """ Computes the wave length at given x,t

    Input:
        data:       dataset containing all data
        t:          t-coordinate
        x:          x-coordinate

    Options:
        crests:     find crests (True) or throughs (default: False)
        no_result:  value to return if there is no result (default: np.nan)

    Output:
        lengths:    list of distances between wave crests
    """
    y_idx = fu.find_peaks_const_t(data, t, x=x, crests=crests)
    y_idx = np.sort(y_idx)

    if y_idx.size < 2:
        print("No wave period can be computed!")
        return np.array([no_result])

    lengths = data["y"][y_idx].diff(dim="y").values

    return lengths


def spectral_analysis_1d(data: xr.Dataset, y: Numeric, x: Numeric = 1e4, variable: str = "wl", demean: bool = False) -> Tuple[np.ndarray]:
    """ Apply fourier transform on timeseries

    Input:
        data:       Dataset containing all data and coordinates
        y:          y-coordinate
        x:          x-coordinate

    Options:
        variable:   name of the variable to use, e.g. "wl" or "p", default="wl"
        demean:     remove mean from signal, default=False

    Output:
        freqs:      corresponding frequencies (time-domain)
        power:      power-spectrum
    """
    t = data["t"].values.astype("datetime64[s]").astype(float)
    dt = np.median(np.diff(t))
    var = data[variable].interp(x=x, y=y)

    if demean:
        var -= np.mean(var)

    transform = np.fft.rfft(var)
    power = np.power(np.abs(transform), 2.)

    freqs = np.fft.rfftfreq(var.size, dt)

    return (freqs, power)


def spectral_analysis_2d(data: xr.Dataset, x: Numeric = 1e4, variable: str = "wl", demean: bool = False) -> Tuple[np.ndarray]:
    """ Apply fourier transform on spatial and temporal varying data

    Input:
        data:       Dataset containing all data and coordinates

    Parameters:
        x:          x-coordinate
        variable:   name of the variable to use, e.g. "wl" or "p", default="wl"
        demean:     remove mean from signal, default=False

    Output:
        wavenumber: corresponding wavenumbers (space-domain)
        freqs:      corresponding frequencies (time-domain)
        power:      power-spectrum
    """
    t = data["t"].values.astype("datetime64[s]").astype(float)
    dt = np.median(np.diff(t))
    y = data["y"].values
    dy = np.median(np.diff(y))
    var = data[variable].interp(x=x)

    if demean:
        var -= np.mean(var)

    # transform first over time, then space
    transform = np.fft.rfft2(var, axes=(1, 0))
    # rearrange data, so that wavenumber is in increasing order
    transform = np.fft.fftshift(transform, axes=1)
    power = np.power(np.abs(transform), 2.)

    freqs = np.fft.rfftfreq(var.shape[0], dt)
    wavenumber = np.fft.fftfreq(var.shape[1], dy)
    wavenumber = np.fft.fftshift(wavenumber)

    return wavenumber, freqs, power
