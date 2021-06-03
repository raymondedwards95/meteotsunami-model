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

    ## Figure options
    fig, ax = plt.subplots(1, 2, sharey=True)
    fig.set_size_inches(14.4, 7.2)
    fig.set_dpi(100)
    fig.set_tight_layout(True)

    ## Initial data
    plotdata = np.zeros(2, dtype=np.object)
    plottext = np.zeros(2, dtype=np.object)

    plotdata[0] = ax[0].contourf(
        dataset["x"] / 1000,
        dataset["y"] / 1000,
        dataset["wl"].isel(t=0)
    )
    plotdata[1] = ax[1].contourf(
        dataset["x"] / 1000,
        dataset["y"] / 1000,
        dataset["p"].isel(t=0),
        vmin=np.floor(dataset["p"].min()),
        vmax=np.ceil(dataset["p"].max())
    )
    plottext[0] = ax[0].set_title(
        f"Sea Surface Elevation  \n$t={dataset['t'].isel(t=0).values.tolist()/1e9/3600}$ hours since start"
    )
    plottext[1] = ax[1].set_title(
        f"Surface Air Pressure  \n$t={dataset['t'].isel(t=0).values.tolist()/1e9/3600}$ hours since start"
    )

    ## Subplot layout
    def initfig():
        for i in range(2):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_xlabel("$x$ [km]")
            ax[i].set_xlim([0, 200])
            ax[i].set_ylim([0, dataset["y"].max()])
        
        ax[0].set_ylabel("$y$ [km]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())
    
    initfig()

    ## Update data
    def update(i):
        plotdata[0].set_array(dataset["wl"].isel(t=i))
        plotdata[1].set_array(dataset["p"].isel(t=i))
        plottext[0].set_text(
            f"Sea Surface Elevation  \n$t={dataset['t'].isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )
        plottext[1].set_text(
            f"Surface Air Pressure  \n$t={dataset['t'].isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())
    
    ## Animation
    frames = (np.arange(dataset["t"].size)).astype(np.int)
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
