""" Functions for visualising model outputs """

import datetime
import os
import sys
import time

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.animation import FuncAnimation


@np.vectorize
def to_timestr(seconds):
    """ Converts time in seconds since reference to a date-string """
    return datetime.datetime.fromtimestamp(seconds).strftime("%Y-%m-%d %H:%M:%S")


def animation_contour(dataset, saveloc=None):
    """ Creates an animation of the top-down view of the water level and surface air pressure data

    Input:
        dataset:    Dataset that contains all variables
    
    Parameters:
        saveloc:    Folder to write the animation to
    """
    assert isinstance(dataset, xr.Dataset)

    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".mp4"):
        saveloc.replace(".mp4", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + "/anim_contours.mp4"
    savename_static = saveloc + "/test_anim_contours.jpg"

    ## Shortcuts
    x = dataset["x"]
    y = dataset["y"]
    t = dataset["t"]
    wl = dataset["wl"]
    p = dataset["p"]

    x = x - x.min()
    wl_max = np.max([np.abs([wl.max(), wl.min()])])
    wl_min = -1. * wl_max
    p_max = np.max([np.abs([p.max(), p.min()])])
    p_min = -1. * p_max

    limits = [(wl_min, wl_max), (p_min, p_max)]
    cmaps = [cmo.cm.balance, cmo.cm.delta]

    ## Figure options
    fig, ax = plt.subplots(2, 1, sharey=True)
    fig.set_size_inches(14.4, 7.2)
    fig.set_dpi(100)
    fig.set_tight_layout(True)

    ## Initial data
    plotdata = np.zeros(2, dtype=np.object)
    plottext = np.zeros(1, dtype=np.object)

    def set_plotdata(i=0):
        plotdata[0] = ax[0].contourf(
            y / 1000,
            x / 1000,
            wl.isel(t=i).T,
            vmin=wl_min,
            vmax=wl_max,
            cmap=cmaps[0],
            levels=31
        )
        plotdata[1] = ax[1].contourf(
            y / 1000,
            x / 1000,
            p.isel(t=i).T,
            vmin=p_min,
            vmax=p_max,
            cmap=cmaps[1],
            levels=31
        )

    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )

    set_plotdata()
    set_plottext()

    ## Subplot layout
    def initfig():
        for i in range(2):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_ylabel("$x$ [km]")
            ax[i].set_ylim([0, 200])
            ax[i].set_xlim([0, y.max() / 1000.])

            # colorbar
            plotdata[i].set_clim(limits[i])
            cbar = fig.colorbar(
                plotdata[i],
                ax=ax[i],
                fraction=0.01,
                aspect=50
            )

        ax[0].set_xlabel("$y$ [km]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    initfig()

    ## Update data
    num_frames = t.size

    def update(i):
        # progress
        if not (num_frames-i-1) % (num_frames // 5):
            print(f"Frame {i:4.0f} of {num_frames:0.0f} ({(i+1)/num_frames*100:0.1f}%)")

        # remove data in contour plots
        for _temp in plotdata[0].collections:
            _temp.remove()
        for _temp in plotdata[1].collections:
            _temp.remove()

        # new data
        set_plotdata(i)
        set_plottext(i)

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    ## Animation
    frames = (np.arange(t.size)).astype(np.int)
    anim = FuncAnimation(
        fig,
        update,
        init_func=initfig,
        frames=frames,
        interval=1000/20
    )
    print(f"Creating animation '{savename}'")
    t0 = time.perf_counter()
    plt.savefig(savename_static)
    anim.save(savename)
    t1 = time.perf_counter()
    print(f"Finished contour-animation in {t1-t0:0.1f} seconds")
    print(f"Saved animation as '{savename}'")
    return


def animation_alongshore(dataset, saveloc=None):
    """ Creates an animation of an alongshore cross-section of the water level and surface air pressure data

    Input:
        dataset:    Dataset that contains all variables
    
    Parameters:
        saveloc:    Folder to write the animation to
    """
    assert isinstance(dataset, xr.Dataset)

    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".mp4"):
        saveloc.replace(".mp4", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + "/anim_alongshore.mp4"
    savename_static = saveloc + "/test_anim_alongshore.jpg"

    ## Shortcuts
    x = dataset["x"]
    y = dataset["y"]
    t = dataset["t"]
    wl = dataset["wl"]
    p = dataset["p"]

    x = x - x.min()
    wl_max = np.max([np.abs([wl.max(), wl.min()])])
    wl_min = -1. * wl_max
    p_max = np.max([np.abs([p.max(), p.min()])])
    p_min = -1. * p_max

    ## Figure options
    fig, ax = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(14.4, 7.2)
    fig.set_dpi(100)
    fig.set_tight_layout(True)

    ## Initial data
    plotdata = np.zeros(2, dtype=np.object)
    plottext = np.zeros(1, dtype=np.object)

    def set_plotdata(i=0):
        plotdata[0] = ax[0].plot(
            y / 1000,
            wl.isel(t=i).interp(x=10000),
            color="C0"
        )[0]
        plotdata[1] = ax[1].plot(
            y / 1000,
            p.isel(t=i).interp(x=10000),
            color="C1"
        )[0]
    
    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )
    
    set_plotdata()
    set_plottext()

    ## Subplot layout
    def initfig():
        for i in range(2):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_xlim([y.min() / 1000. / 10., y.max() / 1000.])
        
        ax[0].set_ylim([wl_min, wl_max])
        ax[1].set_ylim([p_min, p_max])
        ax[1].set_xlabel("$y$ [km]")
        ax[0].set_ylabel("Sea Surface Elevation [m]")
        ax[1].set_ylabel("Surface Air Pressure [Pa]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())
    
    initfig()

    ## Update data
    num_frames = t.size
    def update(i):
        # progress
        if not (num_frames-i-1) % (num_frames // 5):
            print(f"Frame {i:4.0f} of {num_frames:0.0f} ({(i+1)/num_frames*100:0.1f}%)")

        # new data
        plotdata[0].set_ydata(wl.isel(t=i).interp(x=10000))
        plotdata[1].set_ydata(p.isel(t=i).interp(x=10000))
        set_plottext(i)

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())
    
    ## Animation
    frames = (np.arange(t.size)).astype(np.int)
    anim = FuncAnimation(
        fig,
        update,
        init_func=initfig,
        frames=frames,
        interval=1000/20
    )
    print(f"Creating animation '{savename}'")
    t0 = time.perf_counter()
    plt.savefig(savename_static)
    anim.save(savename)
    t1 = time.perf_counter()
    print(f"Finished alongshore-animation in {t1-t0:0.1f} seconds")
    print(f"Saved animation as '{savename}'")
    return


def animation_crossshore(dataset, saveloc=None):
    """ Creates an animation of crossshore cross-sections of the water level and surface air pressure data

    Input:
        dataset:    Dataset that contains all variables
    
    Parameters:
        saveloc:    Folder to write the animation to
    """
    assert isinstance(dataset, xr.Dataset)

    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".mp4"):
        saveloc.replace(".mp4", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + "/anim_crossshore.mp4"
    savename_static = saveloc + "/test_anim_crossshore.jpg"

    ## Shortcuts
    x = dataset["x"]
    y = dataset["y"]
    t = dataset["t"]
    wl = dataset["wl"]

    x = x - x.min()
    wl_max = np.max([np.abs([wl.max(), wl.min()])])
    wl_min = -1. * wl_max
    y_max = y.max().values

    ## Figure options
    slices = 5
    fig, ax = plt.subplots(slices, 1, sharex=True)
    fig.set_size_inches(14.4, 7.2)
    fig.set_dpi(100)
    fig.set_tight_layout(True)

    ## Initial data
    plotdata = np.zeros(slices, dtype=np.object)
    plottext = np.zeros(1, dtype=np.object)
    _y = np.arange(slices) * y_max / slices

    def set_plotdata(i=0):
        for j in range(slices):
            plotdata[j] = ax[j].plot(
                x / 1000,
                wl.isel(t=i).interp(y=_y[j]),
                color="C0",
                label=f"$y={_y[j]/1000:0.0f}$ km"
            )[0]
    
    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )
    
    set_plotdata()
    set_plottext()

    ## Subplot layout
    def initfig():
        for i in range(slices):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_xlim([0, 200])
            ax[i].set_ylabel("SSE [m]")
            ax[i].set_ylim([wl_min, wl_max])
            ax[i].legend()
            
        ax[-1].set_xlabel("$x$ [km]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())
    
    initfig()

    ## Update data
    num_frames = t.size
    def update(i):
        # progress
        if not (num_frames-i-1) % (num_frames // 5):
            print(f"Frame {i:4.0f} of {num_frames:0.0f} ({(i+1)/num_frames*100:0.1f}%)")

        # new data
        for j in range(slices):
            plotdata[j].set_ydata(wl.isel(t=i).interp(y=_y[j]))
        set_plottext(i)

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())
    
    ## Animation
    frames = (np.arange(t.size)).astype(np.int)
    anim = FuncAnimation(
        fig,
        update,
        init_func=initfig,
        frames=frames,
        interval=1000/20
    )
    print(f"Creating animation '{savename}'")
    t0 = time.perf_counter()
    plt.savefig(savename_static)
    anim.save(savename)
    t1 = time.perf_counter()
    print(f"Finished crossshore-animation in {t1-t0:0.1f} seconds")
    print(f"Saved animation as '{savename}'")
    return
