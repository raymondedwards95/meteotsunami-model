""" Additional analysis of waterlevel data """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.utilities as fu


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


def _test_compute_decay_parameter(showfig=True):
    """ Test for computation of decay parameter """
    k0 = 1./12.
    y0 = 0.7

    x = np.arange(100)
    y = exp_decay(x, k0, y0) + np.random.normal(scale=0.01, size=x.size)

    (k0_est, y0_est), _ = curve_fit(
        exp_decay,
        x,
        y
    )

    if showfig:
        plt.figure()
        plt.scatter(x, y)
        plt.plot(x, exp_decay(x, k0_est, y0_est), "C1", label=f"$1/k_0={1/k0_est:0.1f}$; $y_0={y0_est:0.2f}$")
        plt.legend()
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.show()
    return


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


def spectral_analysis_1d(data, y, x=1e4, variable="wl"):
    """ Apply fourier transform on timeseries 
    
    Input:
        data:       Dataset containing all data and coordinates
        y:          y-coordinate
    
    Parameters:
        x:          x-coordinate
        variable:   name of the variable to use, e.g. "wl" or "p"
    
    Output:
        freqs:      corresponding frequencies (time-domain)      
        power:      power-spectrum
    """
    t = data["t"].values.astype("datetime64[s]").astype(float)
    dt = np.median(np.diff(t))
    var = data[variable].interp(x=x, y=y)

    transform = np.fft.rfft(var)
    power = np.power(np.abs(transform), 2.)

    freqs = np.fft.rfftfreq(var.size, dt)

    return freqs, power


def spectral_analysis_2d(data, x=1e4, variable="wl"):
    """ Apply fourier transform on spatial and temporal varying data 
    
    Input:
        data:       Dataset containing all data and coordinates
    
    Parameters:
        x:          x-coordinate
        variable:   name of the variable to use, e.g. "wl" or "p"
    
    Output:
        freqs:      corresponding frequencies (time-domain)      
        wavenumber: corresponding wavenumbers (space-domain)
        power:      power-spectrum
    """
    t = data["t"].values.astype("datetime64[s]").astype(float)
    dt = np.median(np.diff(t))
    y = data["y"].values
    dy = np.median(np.diff(y))
    var = data[variable].interp(x=x)

    transform = np.fft.rfft2(var, axes=(1, 0))  # first over time, then space
    transform = np.fft.fftshift(transform, axes=1)  # rearrange data, so that wavenumber is in increasing order
    power = np.power(np.abs(transform), 2.)

    freqs = np.fft.rfftfreq(var.shape[0], dt)
    wavenumber = np.fft.fftfreq(var.shape[1], dy)
    wavenumber = np.fft.fftshift(wavenumber)

    return wavenumber, freqs, power


if __name__ == "__main__":
    _test_compute_decay_parameter(showfig=True)
