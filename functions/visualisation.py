""" Functions for visualising model outputs """

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

def vis_timeseries(data, x=1e4, y=1e5):
    if isinstance(y, (int, float)):
        y = [y]

    t = data["t"]
    wl = data["wl"].interp(x=x, y=y)
    p = data["p"].interp(x=x, y=y)

    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(t, wl)
    ax[1].plot(t, p)
    for i in range(2):
        ax[i].axhline(color="black", linewidth=1)

    fig, ax = plt.subplots(len(y), 1, sharex=True)
    for i in range(len(y)):
        ax[i].plot(t, wl[:, i])
        ax[i].plot(t, p[:, i] / 2000.)
        ax[i].axhline(color="black", linewidth=1)

    fig, ax = plt.subplots(2, len(y), sharex=True, sharey=True)
    for i in range(len(y)):
        ax[0,i].plot(t, wl[:, i])
        ax[1,i].plot(t, p[:, i] / 2000.)
        for j in range(2):
            ax[j,i].axhline(color="black", linewidth=1)


def vis_alongshore(data, t=3600, x=1e4):
    raise NotImplementedError


def vis_crossshore(data, y=1e5, t=3600):
    raise NotImplementedError
