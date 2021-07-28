""" Functions for visualising model outputs """

import os
import sys

import cmocean as cmo
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.analysis as fa
import functions.utilities as fu


def vis_timeseries(data, y, x=1e4, t_max=None, saveloc=None, keep_open=False):
    ## Check input
    if not np.isscalar(x):
        raise ValueError(f"{x=} should be scalar")

    if t_max is None:
        t_max = np.inf
    elif np.isscalar(t_max):
        t_max /= 3600

    if np.isscalar(y):
        y = np.array([y])
    if isinstance(y, (list, np.ndarray)):
        y = np.array(y)
    
    y_arr = y
    del y

    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/timeseries_{y_arr.size}"

    ## Extract data
    t = data["t"].values.astype("datetime64[s]").astype(float) / 3600
    wl = data["wl"].interp(x=x, y=y_arr)
    p = data["p"].interp(x=x, y=y_arr)

    t_max = np.min([t.max(), t_max])
    y_max = np.max(data["y"].values)
    wl_max = np.max(np.abs(wl.values))
    p_max = np.max(np.abs(p.values))

    y_arr = np.array([y for y in y_arr if y < y_max])
    if y_arr.size < 1:
        raise ValueError("Input 'y' has no values on the domain of data")
    
    wl_t_idx = [fu.find_peaks_const_y(data, y, x=x, crests=False, variable="wl") for y in y_arr]
    p_t_idx = [fu.find_peaks_const_y(data, y, x=x, crests=True, variable="p") for y in y_arr]

    ## Figure a
    fig, ax = plt.subplots(2, 1, sharex=True, squeeze=False)
    fig.set_size_inches(14, 7)
    fig.set_dpi(FIG_DPI)
    # fig.set_tight_layout(True)

    ax = np.ravel(ax)
    ax[0].plot(t, wl)
    ax[1].plot(t, p)

    for i in range(y_arr.size):
        try:
            for j in wl_t_idx[i]:
                ax[0].axvline(t[j], color=f"C{i}", alpha=0.5)
                ax[0].annotate(f"$t={t[j]:0.1f}$h", xy=(t[j], wl_max), ha="left", va="top", color=f"C{i}", rotation=90)
        except:
            pass
        try:
            for j in p_t_idx[i]:
                ax[1].axvline(t[j], color=f"C{i}", alpha=0.5)
                ax[1].annotate(f"$t={t[j]:0.1f}$h", xy=(t[j], 0), ha="left", va="bottom", color=f"C{i}", rotation=90)
        except:
            pass

    for i in range(2):
        ax[i].axhline(color="black", linewidth=1)
        ax[i].set_xlim([0, t_max])
        ax[i].grid()
        ax[i].xaxis.set_minor_locator(MultipleLocator(1))

    ax[-1].set_xlabel("Time since start [hours]")
    ax[0].set_ylabel("Sea Surface Elevation [m]")
    ax[1].set_ylabel("Atmospheric Pressure [Pa]")
    ax[0].legend([f"$y={y/1000:0.0f}$km" for y in y_arr], bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=4, mode="expand", borderaxespad=0.)

    fig.savefig(savename + "_a", bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure {savename}_a")

    ## Figure b
    fig, ax = plt.subplots(y_arr.size, 1, sharex=True, squeeze=False)
    fig.set_size_inches(FIGSIZE_LONG)
    fig.set_dpi(FIG_DPI)
    ax = np.ravel(ax)
    ax2 = np.array([ax[i].twinx() for i in range(ax.size)])
    for i in range(y_arr.size):
        wl_slice = wl[:, i].values
        p_slice = p[:, i].values
        ax[i].plot(t, wl_slice, color="C0")
        ax2[i].plot(t, p_slice, color="C1")
        ax[i].axhline(color="black", linewidth=1)

        ax[i].set_xlim([0, t_max])
        ax[i].set_ylim(np.array([-1.1, 1.1]) * np.max(np.abs(data["wl"].values)))
        ax2[i].set_ylim(np.array([-1.1, 1.1]) * np.max(np.abs(data["p"].values)))
        ax[i].xaxis.set_minor_locator(MultipleLocator(1))

        ax[i].tick_params(axis="y", labelcolor=f"C0")
        ax2[i].tick_params(axis="y", labelcolor=f"C1")

        ax[i].grid()
        ax[i].annotate(f"$y={y_arr[i]/1000:0.0f}$km", xy=(0.01, 0.99), xycoords="axes fraction", ha="left", va="top")

        try:
            for j in wl_t_idx[i]:
                ax[i].axvline(t[j], color="C0", alpha=0.5)
                ax[i].annotate(f"$t={t[j]:0.1f}$h", xy=(t[j], wl_max), ha="left", va="top", color="C0", rotation=90)
        except:
            pass
        try:
            for j in p_t_idx[i]:
                ax[i].axvline(t[j], color="C1", alpha=0.5)
                ax[i].annotate(f"$t={t[j]:0.1f}$h", xy=(t[j], 0), ha="left", va="bottom", color="C1", rotation=90)
        except:
            pass

    ax[-1].set_ylabel("Water Level [m]")
    ax2[-1].set_ylabel("Surface Pressure [Pa]")
    ax[-1].set_xlabel("Time since start [hours]")

    fig.savefig(savename + "_b", bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure {savename}_b")

    ## End
    if not keep_open:
        plt.close("all")


def vis_alongshore(data, t=3600, x=1e4, saveloc=None, keep_open=False):
    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/along_shore_{x/1000:0.0f}_{t/3600:0.0f}"

    ## Check input
    # if np.isscalar(t):
    #     t = np.array([t])
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")

    # assert t.size == 1

    ## Extract data
    y = data["y"]
    wl = data["wl"].interp(t=fu.to_timestr(t), x=x)

    ## Figure
    fig, ax = plt.subplots(1, 1, squeeze=False)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    ax = np.ravel(ax)

    ax[0].plot(
        y / 1000., 
        wl
    )
    ax[0].axhline(color="black", linewidth=1)
    ax[0].set_xlim(0, y.max() / 1000.)
    ax[0].set_xlabel("$y$ [km]")
    ax[0].set_ylabel("$SSE$ [m]")
    ax[0].set_title(f"Along-shore profile at $x={x/1000:0.0f}$km and $t={t/3600:0.1f}$hours")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")


def vis_crossshore(data, y=1e5, t=3600, saveloc=None, keep_open=False):
    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/cross_shore_{y/1000:0.0f}_{t/3600:0.0f}"

    ## Check input
    # if np.isscalar(y):
    #     y = np.array([y])
    # if np.isscalar(t):
    #     t = np.array([t])

    # assert y.size == 1
    # assert t.size == 1

    ## Make best fit
    k0, y0 = fa.compute_decay_parameter(data, y, t)
    wl_model = fa.exp_decay(data["x"], k0, y0)

    ## Extract data
    x = data["x"]
    wl = data["wl"].interp(y=y, t=fu.to_timestr(t))

    ## Figure
    fig, ax = plt.subplots(1, 1, squeeze=False)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    ax = np.ravel(ax)
    ax[0].plot(
        x / 1000., 
        wl, 
        color="C0", 
        label="Waterlevel"
    )
    ax[0].plot(
        x / 1000., 
        wl_model, 
        color="C0", 
        linestyle="--", 
        label=f"Best fit with $1/k0={1./k0/1000.:0.1f}$km"
    )
    ax[0].axhline(color="black", linewidth=1)
    ax[0].legend()
    ax[0].set_xlim(0, x.max() / 1000.)
    ax[0].set_xlabel("$x$ [km]")
    ax[0].set_ylabel("$SSE$ [m]")
    ax[0].set_title(f"Cross-shore profile at $y={y/1000:0.0f}$km and $t={t/3600:0.1f}$hours")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")


def vis_spectrum_1d(data, x=1e4, y=1e5, saveloc=None, keep_open=False, variable="wl", demean=True):
    ## Check inputs
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")

    if np.isscalar(y):
        figsize = FIGSIZE_NORMAL
        y_list = np.array([y])
    elif isinstance(y, (list, np.ndarray)):
        figsize = FIGSIZE_HIGH
        y_list = np.array(y)
    else:
        raise ValueError(f"{y=} is not a number or list of numbers")

    variable = variable.lower()
    if variable == "wl":
        variable_name = "Water Level"
    elif variable == "u":
        variable_name = "x-component of velocity"
    elif variable == "v":
        variable_name = "y-component of velocity"
    elif variable == "p":
        variable_name = "Surface Air Pressure"
    else:
        raise ValueError(f"{variable=} should be 'wl', 'u', 'v' or 'p'")

    del y
    y_list = np.ravel(y_list)
    y_num = y_list.size

    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    if y_num == 1:
        savename = saveloc + f"/spectrum_1d_{x/1000:0.0f}_{y_list[0]/1000:0.0f}_{variable}_{int(demean)}"
    else:
        savename = saveloc + f"/spectrum_1d_{x/1000:0.0f}_n{y_num:0.0f}_{variable}_{int(demean)}"

    ## Compute spectrum
    freqs_all = [None] * y_num
    power_all = [None] * y_num
    for i in range(y_num):
        freqs_all[i], power_all[i] = fa.spectral_analysis_1d(data, y_list[i], x=x, variable=variable, demean=demean)

    ## Figure
    fig, ax = plt.subplots(y_num, 1, squeeze=False, sharex=True)
    fig.set_size_inches(figsize)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    fig.suptitle(f"Power Spectrum - {variable_name}", va="top", ha="left", x=0.01)
    ax = np.ravel(ax)

    for i in range(y_num):
        ax[i].plot(
            freqs_all[i] * 3600.,
            power_all[i] / 3600.,
            color=f"C{i}"
        )
        ax[i].fill_between(
            freqs_all[i] * 3600, 
            power_all[i] / 3600, 
            alpha=0.1,
            color=f"C{i}"
        )
        ax[i].axhline(color="black", linewidth=1)
        ax[i].set_ylim(0, None)
        ax[i].set_xlim(0, 2.2)
        ax[i].set_ylabel("Spectral Power [m$^2$ hr]")
        ax[i].legend([f"$y = {y_list[i]/1000:0.0f}$km"], loc="upper right")
        ax[i].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax[i].xaxis.set_minor_locator(MultipleLocator(0.1))
        ax[i].grid()
    ax[-1].set_xlabel("Frequency [cycles / hour]")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")
    

def vis_spectrum_2d(data, x=1e4, saveloc=None, keep_open=False, variable="wl", xlims=None, ylims=None, autolim=1e-3):
    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/spectrum_2d_{x/1000:0.0f}_{variable}"

    ## Check inputs
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")

    ## Compute spectrum
    wavenumber, freqs, power = fa.spectral_analysis_2d(data, x=x, variable=variable)

    ## Compute plot limits
    if (xlims is None) or (ylims is None):
        valid = (power / power.max()) > autolim
    if xlims is None:
        i = np.any(valid, axis=0)
        xmax = 1e6 * np.min([wavenumber.max(), 1.5*np.max(np.abs(wavenumber[i]))])
        xlims = [0, xmax]
    if ylims is None:
        i = np.any(valid, axis=1)
        ymax = 3600. * np.min([freqs.max(), 1.5*np.max(freqs[i])])
        ylims = [0, ymax]

    ## Figure
    fig, ax = plt.subplots(1, 1, squeeze=False)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    ax = np.ravel(ax)
    div = make_axes_locatable(ax[0])
    cax = div.append_axes("right", "5%", "5%")

    im = ax[0].pcolormesh(  # note: contourf is an option
        np.flip(wavenumber) * 1e6,
        freqs * 3600.,
        power,
        shading="nearest",
        cmap=cmo.cm.matter,
        vmin=5.*power.min()
    )
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Spectral Power")

    for i in range(5):
        ax[0].plot(
            wavenumber * 1e6, 
            fa.dispersion_relation(wavenumber, n=i, alpha=1/400) * 3600., 
            linewidth=1, 
            color="black"
        )
    ax[0].axvline(color="black", linewidth=1)
    ax[0].set_xlim(xlims)
    ax[0].set_ylim(ylims)
    ax[0].set_xlabel("Wavenumber [1 / 1000km]")
    ax[0].set_ylabel("Frequency [cycles / hour]")
    ax[0].set_title(f"Power Spectrum at $x={x/1000:0.0f}$km")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")
    

def vis_contour(data, t, saveloc=None, keep_open=False, variable="wl", xlims=None, ylims=None):
    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/contour_{variable}"

    ## Check input
    if variable in ["wl"]:
        cmap = cmo.cm.balance
    elif variable in ["u", "v"]:
        cmap = cmo.cm.delta
    else:
        raise ValueError(f"Excpected {variable=} to be 'wl', 'u' or 'v'")

    if np.isscalar(t):
        t_arr = np.array([t])
    elif isinstance(t, list):
        t_arr = np.array(t)
    elif isinstance(t, (list, np.ndarray)):
        t_arr = t
    else:
        raise ValueError(f"{t=} is not an array, but {type(t)}")
    
    t_num = t_arr.size
    del t

    ## Get parameters
    x = data["x"]
    y = data["y"]
    var = data[variable]

    x = x - x.min()
    var_max = float(np.max([np.abs([var.max(), var.min()])]))
    var_min = -1. * var_max

    if xlims is None:
        xlims = [y.min() / 1000. / 10., y.max() / 1000.]

    if ylims is None:
        ylims = [0, 200]

    ## Figure
    fig, ax = plt.subplots(t_num, 1, squeeze=False, sharex=True)
    fig.set_size_inches(FIGSIZE_LONG)
    fig.set_dpi(FIG_DPI)
    # fig.set_tight_layout(True)
    ax = np.ravel(ax)

    norm = Normalize(var_min, var_max)
    im = cm.ScalarMappable(norm=norm, cmap=cmap)

    ## Subplots
    for i in range(t_num):
        ax[i].contourf(
            y / 1000, 
            x / 1000, 
            var.interp(t=fu.to_timestr(t_arr[i])).T, 
            100,
            cmap=cmap,
            norm=norm,
            vmin=var_min,
            vmax=var_max    
        )
        ax[i].set_ylabel("$x$ [km]")
        # ax[i].set_title(f"$t = {t_arr[i] / 3600:0.1f}$h")
        ax[i].annotate(
            f"$t = {t_arr[i] / 3600:0.1f}$h", 
            xy=(0.99, 0.98),
            xycoords="axes fraction",
            ha="right", 
            va="top"
        )
        ax[i].set_xlim(xlims)
        ax[i].set_ylim(ylims)
    
    ax[-1].set_xlabel("$y$ [km]")
    fig.colorbar(
        im, 
        ax=ax.ravel().tolist(),
        label="Sea Surface Elevation [m]",
        aspect=50
    )

    ## Save figure
    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")


