""" Functions for visualising model outputs """

import datetime
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.animation import FuncAnimation


@np.vectorize
def to_timestr(seconds):
    """ Converts time in seconds since reference to a date-string """
    return datetime.datetime.fromtimestamp(seconds).strftime("%Y-%m-%d %H:%M:%S")


def animation_contour(dataset, saveloc=None):
    assert isinstance(dataset, xr.Dataset)

    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".mp4"):
        saveloc.replace(".mp4", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + "/anim_contours.mp4"

    ## Shortcuts
    x = dataset["x"]
    y = dataset["y"]
    t = dataset["t"]
    wl = dataset["wl"]
    p = dataset["p"]

    x = x - x.min()
    wl_max = np.max([wl.max(), np.abs(wl.min())])
    p_max = np.ceil(p.max())
    p_min = np.floor(p.min())

    ## Figure options
    fig, ax = plt.subplots(1, 2, sharey=True)
    fig.set_size_inches(14.4, 7.2)
    fig.set_dpi(100)
    fig.set_tight_layout(True)

    ## Initial data
    plotdata = np.zeros(2, dtype=np.object)
    plottext = np.zeros(2, dtype=np.object)

    def set_plotdata(i=0):
        plotdata[0] = ax[0].contourf(
            x / 1000,
            y / 1000,
            wl.isel(t=i),
            vmin=-1.*wl_max,
            vmax=wl_max
        )
        plotdata[1] = ax[1].contourf(
            x / 1000,
            y / 1000,
            p.isel(t=i),
            vmin=p_min,
            vmax=p_max
        )
    
    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"Sea Surface Elevation  \n$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )
        plottext[1] = ax[1].set_title(
            f"Surface Air Pressure  \n$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )
    
    set_plotdata()
    set_plottext()

    ## Subplot layout
    def initfig():
        for i in range(2):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_xlabel("$x$ [km]")
            ax[i].set_xlim([0, 200])
            ax[i].set_ylim([0, y.max() / 1000.])
        
        ax[0].set_ylabel("$y$ [km]")
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
    anim.save(savename)
    t1 = time.perf_counter()
    print(f"Finished contour-animation in {t1-t0:0.1f} seconds")
    print(f"Saved animation as '{savename}'")
    return
