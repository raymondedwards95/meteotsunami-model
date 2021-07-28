""" Additional analysis of waterlevel data """

import os
import sys

import numpy as np
from scipy.optimize import curve_fit

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.utilities as fu


def dispersion_relation(k, n=0, alpha=1/400):
    """ Returns the frequency related to the wavenumber 
    for shallow water waves on a beach of linear slope [Eckart, 1951] 
    
    Input:
        k:      wavenumbers
    
    Parameters:
        n:      mode
        alpha:  slope of beach
    
    Output:
        f:      frequencies corresponding to wavenumbers
    """
    g = 9.81
    return np.sqrt((g * np.abs(k) / 2. / np.pi) * (2. * n + 1.) * np.tan(alpha))


def exp_decay(x, k0, y0):
    """ Returns an exponential decay in x with decay parameter k0 and scale y0 """
    return y0 * np.exp(-k0 * x)


def compute_decay_parameter(data, y, t):
    """"Computes the decay parameter and scaling parameter for data at given y and t """
    popt, _ = curve_fit(
        exp_decay,  # f
        data["x"],  # xdata
        data["wl"].interp(t=fu.to_timestr(t), y=y),  # ydata
        p0=[1./100000., 1.],  # p0
    )
    k0, y0 = popt
    return k0, y0


def compute_wave_periods(data, y, x=None, crests=True):
    """ Computes the wave period at given x,y 
    
    Input:
        data:       dataset containing all data
        y:          y-coordinate
    
    Parameters:
        x:          x-coordinate
        crests:     find crests (True) or throughs (False)
    
    Output:
        periods:    list of time-intervals between wave crests
    """
    t_idx = fu.find_peaks_const_y(data, y, x=x, crests=crests)
    t_idx = np.sort(t_idx)

    if t_idx.size < 2:
        print("No wave period can be computed!")
        return np.array([np.nan])

    periods = data["t"][t_idx].diff(dim="t").values.astype("timedelta64[s]").astype(float) / 3600

    return periods


def compute_wave_lengths(data, t, x=None, crests=True):
    """ Computes the wave length at given x,t 
    
    Input:
        data:       dataset containing all data
        t:          t-coordinate
    
    Parameters:
        x:          x-coordinate
        crests:     find crests (True) or throughs (False)
    
    Output:
        lengths:    list of distances between wave crests
    """
    y_idx = fu.find_peaks_const_t(data, t, x=x, crests=crests)
    y_idx = np.sort(y_idx)

    if y_idx.size < 2:
        print("No wave period can be computed!")
        return np.array([np.nan])

    lengths = data["y"][y_idx].diff(dim="y").values

    return lengths


def spectral_analysis_1d(data, y, x=1e4, variable="wl", demean=False):
    """ Apply fourier transform on timeseries 
    
    Input:
        data:       Dataset containing all data and coordinates
        y:          y-coordinate
    
    Parameters:
        x:          x-coordinate
        variable:   name of the variable to use, e.g. "wl" or "p"
        demean:     remove mean from signal
    
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

    return freqs, power


def spectral_analysis_2d(data, x=1e4, variable="wl", demean=False):
    """ Apply fourier transform on spatial and temporal varying data 
    
    Input:
        data:       Dataset containing all data and coordinates
    
    Parameters:
        x:          x-coordinate
        variable:   name of the variable to use, e.g. "wl" or "p"
        demean:     remove mean from signal
    
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

    transform = np.fft.rfft2(var, axes=(1, 0))  # first over time, then space
    transform = np.fft.fftshift(transform, axes=1)  # rearrange data, so that wavenumber is in increasing order
    power = np.power(np.abs(transform), 2.)

    freqs = np.fft.rfftfreq(var.shape[0], dt)
    wavenumber = np.fft.fftfreq(var.shape[1], dy)
    wavenumber = np.fft.fftshift(wavenumber)

    return wavenumber, freqs, power
