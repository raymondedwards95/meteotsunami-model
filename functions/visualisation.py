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
    y_num = y_arr.size
    del y

    if y_num == 1:
        figsize_a = (14, 7)
        figsize_b = FIGSIZE_WIDE
    else:
        figsize_a = (14, 7)
        figsize_b = FIGSIZE_LONG

    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/timeseries_x{x/1000:0.0f}_yn{y_num}"

    ## Extract data
    t = data["t"].values.astype("datetime64[s]").astype(float) / 3600
    wl = data["wl"].interp(x=x, y=y_arr)
    p = data["p"].interp(x=x, y=y_arr)

    t_max = np.nanmin([t.max(), t_max])
    y_max = np.nanmax(data["y"].values)
    wl_max = np.nanmax(np.abs(wl.values))
    p_max = np.nanmax(np.abs(p.values))

    y_arr = np.array([y for y in y_arr if y < y_max])
    if y_num < 1:
        raise ValueError("Input 'y' has no values on the domain of data")
    
    wl_t_idx = [fu.find_peaks_const_y(data, y, x=x, crests=False, variable="wl") for y in y_arr]
    p_t_idx = [fu.find_peaks_const_y(data, y, x=x, crests=True, variable="p") for y in y_arr]

    ## Figure a
    fig, ax = plt.subplots(2, 1, sharex=True, squeeze=False)
    fig.set_size_inches(figsize_a)
    fig.set_dpi(FIG_DPI)
    # fig.set_tight_layout(True)

    ax = np.ravel(ax)
    ax[0].plot(t, wl)
    ax[1].plot(t, p)

    for i in range(y_num):
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
    ax[0].set_ylim(np.array([-1.1, 1.1]) * wl_max)
    ax[0].legend([f"$y={y/1000:0.0f}$ km" for y in y_arr], bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=4, mode="expand", borderaxespad=0.)

    fig.savefig(savename + "_a", bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    print(f"Saved figure {savename}_a")

    ## Figure b
    fig, ax = plt.subplots(y_num, 1, sharex=True, squeeze=False)
    fig.set_size_inches(figsize_b)
    fig.set_dpi(FIG_DPI)
    ax = np.ravel(ax)
    ax2 = np.array([ax[i].twinx() for i in range(ax.size)])
    for i in range(y_num):
        wl_slice = wl[:, i].values
        p_slice = p[:, i].values
        ax[i].plot(t, wl_slice, color="C0")
        ax2[i].plot(t, p_slice, color="C1")
        ax[i].axhline(color="black", linewidth=1)

        ax[i].set_xlim([0, t_max])
        ax[i].set_ylim(np.array([-1.1, 1.1]) * wl_max)
        ax2[i].set_ylim(np.array([-1.1, 1.1]) * p_max)
        ax[i].xaxis.set_minor_locator(MultipleLocator(1))
        ax[i].xaxis.set_tick_params(which="both", top=True)

        ax[i].tick_params(axis="y", labelcolor="C0")
        ax2[i].tick_params(axis="y", labelcolor="C1")

        ax[i].set_ylabel("Water Level [m]", color="C0")
        ax2[i].set_ylabel("Surface Pressure [Pa]", color="C1")

        ax[i].grid()
        ax[i].annotate(f"$y={y_arr[i]/1000:0.0f}$ km", xy=(0.01, 0.99), xycoords="axes fraction", ha="left", va="top")

        try:
            for j in wl_t_idx[i]:
                ax[i].axvline(t[j], color="C0", alpha=0.5)
                ax[i].annotate(f"$t={t[j]:0.1f}$h", xy=(t[j], 0), ha="left", va="bottom", color="C0", rotation=90)
        except:
            pass
        try:
            for j in p_t_idx[i]:
                ax[i].axvline(t[j], color="C1", alpha=0.5)
                ax[i].annotate(f"$t={t[j]:0.1f}$h", xy=(t[j], 0), ha="right", va="top", color="C1", rotation=90)
        except:
            pass
    ax[-1].set_xlabel("Time since start [hours]")

    fig.savefig(savename + "_b", bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    print(f"Saved figure {savename}_b")

    ## End
    if not keep_open:
        plt.close("all")


def vis_alongshore(data, t=3600, x=1e4, saveloc=None, keep_open=False):
    ## Check input
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")

    if np.isscalar(t):
        t_list = np.array([t])
    elif isinstance(t, (list, np.array)):
        t_list = np.array(t)
    else:
        raise ValueError(f"{t=} is not a number or array_like")

    t_num = t_list.size
    del t

    if t_num < 3:
        figsize = FIGSIZE_NORMAL
    else:
        figsize = FIGSIZE_HIGH

    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/along_shore_x{x/1000:0.0f}_tn{t_num}"

    ## Extract data
    y = data["y"]
    wl = data["wl"].interp(t=fu.to_timestr(t_list), x=x)

    ## Figure
    fig, ax = plt.subplots(t_num, 1, squeeze=False, sharex=True)
    fig.set_size_inches(figsize)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    ax = np.ravel(ax)

    for i in range(t_num):
        ax[i].plot(
            y / 1000., 
            wl[i]
        )
        ax[i].axhline(color="black", linewidth=1)
        ax[i].set_xlim(0, y.max() / 1000.)
        ax[i].set_ylabel("$SSE$ [m]")
        ax[i].set_title(f"$x={x/1000:0.0f}$ km and $t={t_list[i]/3600:0.1f}$ hours")
        ax[i].grid()
    ax[-1].set_xlabel("$y$ [km]")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")


def vis_crossshore(data, y=1e5, t=3600, saveloc=None, keep_open=False):
    ## Check input
    if np.isscalar(y):
        y_list = np.array([y])
    elif isinstance(y, (list, np.ndarray)):
        y_list = np.array(y)
    else:
        raise ValueError(f"{y=} is not a number or array_like")

    if np.isscalar(t):
        t_list = np.array([t])
    elif isinstance(t, (list, np.ndarray)):
        t_list = np.array(t)
    else:
        raise ValueError(f"{t=} is not a number or array_like")
    
    del t, y

    y_num = y_list.size
    t_num = t_list.size

    if (t_num == 1) and (y_num == 1):
        figsize = FIGSIZE_NORMAL
    elif (t_num == 1) and (y_num > 1):
        figsize = FIGSIZE_WIDE
    elif (t_num > 1) and (y_num == 1):
        figsize = FIGSIZE_HIGH
    elif (t_num > 1) and (y_num > 1):
        figsize = FIGSIZE_SQUARE
    else:
        raise ValueError(f"Error in determining figsize for {y_num=} and {t_num=}")

    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)

    savename = saveloc + "/cross_shore_"

    if y_num == 1:
        savename += f"y{y_list[0]/1000:0.0f}_"
    else:
        savename += f"yn{y_num}_"
    
    if t_num == 1:
        savename += f"t{t_list[0]/3600:0.0f}"
    else:
        savename += f"tn{t_num}"
    
    ## Figure
    fig, ax = plt.subplots(t_num, y_num, squeeze=False, sharex=True, sharey=True)
    fig.set_size_inches(figsize)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)

    for i in range(t_num):
        for j in range(y_num):
            # Make best fit
            k0, y0 = fa.compute_decay_parameter(data, y_list[j], t_list[i])
            wl_model = fa.exp_decay(data["x"], k0, y0)

            # Plot data
            ax[i,j].plot(
                data["x"] / 1000., 
                data["wl"].interp(y=y_list[j], t=fu.to_timestr(t_list[i])),
                color="C0", 
                label="Waterlevel"
            )

            # Plot best fit
            ax[i,j].plot(
                data["x"] / 1000.,
                wl_model, 
                color="C0", 
                linestyle="--", 
                label=f"Best fit with $1/k0={1./k0/1000.:0.1f}$ km"
            )

            # 
            ax[i,j].axhline(color="black", linewidth=1)
            ax[i,j].legend()
            ax[i,j].set_xlim(0, data["x"].max() / 1000.)
            ax[i,j].set_title(f"$y={y_list[j]/1000:0.0f}$ km; $t={t_list[i]/3600:0.1f}$ hours")
            ax[i,j].grid()

    for i in range(t_num):
        ax[i,0].set_ylabel("$SSE$ [m]")

    for j in range(y_num):
        ax[-1,j].set_xlabel("$x$ [km]")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
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
        ax[i].legend([f"$y = {y_list[i]/1000:0.0f}$ km"], loc="upper right")
        ax[i].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax[i].xaxis.set_minor_locator(MultipleLocator(0.1))
        ax[i].grid()
    ax[-1].set_xlabel("Frequency [cycles / hour]")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")
    

def vis_spectrum_2d(data, x=1e4, saveloc=None, keep_open=False, variable="wl", xlims=None, ylims=None, autolim=1e-3, demean=True):
    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/spectrum_2d_{x/1000:0.0f}_{variable}_{int(demean)}"

    ## Check inputs
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")

    ## Compute spectrum
    wavenumber, freqs, power = fa.spectral_analysis_2d(data, x=x, variable=variable, demean=demean)

    ## Compute plot limits
    if (xlims is None) or (ylims is None):
        valid = (power / power.max()) > autolim
    if xlims is None:
        i = np.any(valid, axis=0)
        xmax = 1e6 * np.nanmin([wavenumber.max(), 1.5*np.nanmax(np.abs(wavenumber[i]))])
        xlims = [0, xmax]
    if ylims is None:
        i = np.any(valid, axis=1)
        ymax = 3600. * np.nanmin([freqs.max(), 1.5*np.nanmax(freqs[i])])
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
    ax[0].set_title(f"Power Spectrum at $x={x/1000:0.0f}$ km")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")
    

def vis_contour(data, t, saveloc=None, keep_open=False, variable="wl", xlims=None, ylims=None):
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
    
    ## Paths
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".jpg"):
        saveloc.replace(".jpg", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/contour_{variable}_n{t_num}_{np.mean(t_arr)/3600:0.0f}"

    ## Get parameters
    x = data["x"]
    y = data["y"]
    var = data[variable]

    x = x - x.min()
    var_max = float(np.nanmax([np.abs([var.max(), var.min()])]))
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
    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    print(f"Saved figure {savename}")

    ## End
    if not keep_open:
        plt.close("all")


if __name__ == "__main__":
    mainpath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    data = xr.open_dataset(f"{mainpath}/reproduction-an-2012/output/data_repr_17.nc")

    vis_contour(data, t=5*3600, keep_open=True)
    vis_contour(data, t=[i * 5 * 3600 for i in range(1, 5)], keep_open=True)

    vis_spectrum_2d(data, x=1e4, keep_open=True)

    vis_spectrum_1d(data, x=1e4, y=1e5, keep_open=True)
    vis_spectrum_1d(data, x=1e4, y=[1e5, 3e5], keep_open=True)
    vis_spectrum_1d(data, x=1e4, y=[1e5, 2e5, 3e5], keep_open=True)
    vis_spectrum_1d(data, x=1e4, y=[1e5, 2e5, 3e5, 4e5], keep_open=True)

    vis_timeseries(data, x=1e4, y=1e5, keep_open=True)
    vis_timeseries(data, x=1e4, y=[i * 10e5 for i in range(1, 5)], keep_open=True, t_max=36*3600)

    vis_crossshore(data, y=1e5, t=3600, keep_open=True)
    vis_crossshore(data, y=[1e5, 2e5], t=3600, keep_open=True)
    vis_crossshore(data, y=1e5, t=[3600, 7200], keep_open=True)
    vis_crossshore(data, y=[1e5, 2e5], t=[3600, 7200], keep_open=True)

    vis_alongshore(data, t=3600, x=1e4, keep_open=True)
    vis_alongshore(data, t=[3000, 7200], x=1e4, keep_open=True)
    vis_alongshore(data, t=[3600, 7200, 10800], x=1e4, keep_open=True)
    
    plt.close("all")

