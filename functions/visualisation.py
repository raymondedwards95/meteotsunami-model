""" Functions for visualising model outputs """

import os
import sys
import time

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.analysis as fa
import functions.utilities as fu


def vis_timeseries(data, x=1e4, y=1e5):
    if np.isscalar(y):
        y = np.array([y])
    if isinstance(y, list):
        y = np.array(y)

    t = data["t"]
    wl = data["wl"].interp(x=x, y=y)
    p = data["p"].interp(x=x, y=y)

    fig, ax = plt.subplots(2, 1, sharex=True, squeeze=False)
    ax = np.ravel(ax)
    ax[0].plot(t, wl)
    ax[1].plot(t, p)
    for i in range(2):
        ax[i].axhline(color="black", linewidth=1)

    fig, ax = plt.subplots(y.size, 1, sharex=True, squeeze=False)
    ax = np.ravel(ax)
    for i in range(y.size):
        ax[i].plot(t, wl[:, i])
        ax[i].plot(t, p[:, i] / 2000.)
        ax[i].axhline(color="black", linewidth=1)

    fig, ax = plt.subplots(2, y.size, sharex=True, sharey=True, squeeze=False)
    for i in range(y.size):
        ax[0,i].plot(t, wl[:, i])
        ax[1,i].plot(t, p[:, i] / 2000.)
        for j in range(2):
            ax[j,i].axhline(color="black", linewidth=1)


def vis_alongshore(data, t=3600, x=1e4, saveloc=None):
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".mp4"):
        saveloc.replace(".mp4", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/along_shore_{x/1000:0.0f}_{t/3600:0.0f}"

    # if np.isscalar(t):
    #     t = np.array([t])
    if not np.isscalar(x):
        raise ValueError(f"{x=} is not a number")

    # assert t.size == 1

    y = data["y"]
    wl = data["wl"].interp(t=fu.to_timestr(t), x=x)

    fig, ax = plt.subplots(1, 1, squeeze=False)
    ax = np.ravel(ax)

    ax[0].plot(y, wl)

    fig.savefig(savename, bbox_inches="tight")


def vis_crossshore(data, y=1e5, t=3600, saveloc=None):
    if saveloc is None:
        saveloc = os.path.dirname(os.path.realpath(__file__)) + "/tests"
    if saveloc.endswith(".mp4"):
        saveloc.replace(".mp4", "")
    os.makedirs(saveloc, exist_ok=True)
    savename = saveloc + f"/cross_shore_{y/1000:0.0f}_{t/3600:0.0f}"

    # if np.isscalar(y):
    #     y = np.array([y])
    # if np.isscalar(t):
    #     t = np.array([t])

    # assert y.size == 1
    # assert t.size == 1

    ## Fit
    k0, y0 = fa.compute_decay_parameter(data, y, t)
    wl_model = fa.exp_decay(data["x"], k0, y0)

    ## Data
    x = data["x"]
    wl = data["wl"].interp(y=y, t=fu.to_timestr(t))

    ## Figure
    fig, ax = plt.subplots(1, 1, squeeze=False)
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

    fig.savefig(savename, bbox_inches="tight")
