""" Script to compare different simulations """

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import os
import sys

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.visualisation as fv


### Cases
cases = [
    [0, 1, 2],
    [0, 3, 4, 5],
    [0, 10, 11, 12],
    [0, 15, 16, 17],
    [0, 20, 21],
]

titles = [
    "Maximum Computational Time Step",
    "Implicitness of Numerical Scheme",
    "Spatial Resolution - Pressure",
    "Temporal Resolution - Pressure",
    "Spatial Resolution - Model"
]

num_comparisons = len(cases)
assert len(titles) == num_comparisons


### Paths
file_dir = os.path.dirname(os.path.realpath(__file__))
output_dir = file_dir + "/output/"
figure_dir = file_dir + "/figures/"


### Defining Visualisations
def vis_alongshore(data_list, title, cases, savename):
    ## Parameters
    _ylims = np.array([[0, 4], [1, 6], [2, 8], [3, 10]]) * 1e3
    _times = [4e4, 8e4, 12e4, 16e4]

    ## Figure
    fig, ax = plt.subplots(2, 2, sharex=False, sharey=True)
    fig.set_size_inches(8, 8)
    fig.set_dpi(150)
    fig.set_tight_layout(True)
    fig.suptitle(f"{title}\nAlong-shore Profile of Sea Surface Elevation")

    # Subplots
    for i in range(4):
        _ax = ax[i//2, i%2]  # select subplot
        _ax.set_xlim(_ylims[i])
        _ax.set_ylim([-0.8, 0.8])
        _ax.set_title(f"$t = {_times[i].:0.0f}$s")
        if i//2: _ax.set_xlabel("$y$ [km]")
        if not i%2: _ax.set_ylabel("$SSE$ [m]")

        for j in range(len(data_list)):
            data = data_list[j]
            _ax.plot(
                data["y"]/1000.,
                data["wl"].interp(t=fv.to_timestr(_times[i]), x=10e3),
                label=f"Case {cases[j]:02.0f}"
            )
        
        if i == 3:
            _ax.legend()
        
    fig.savefig(savename + "/along.jpg", bbox_inches="tight")
    return


def vis_crossshore(data_list, title, cases, savename):
    raise NotImplementedError


def make_comparison(cases, title, id="test"):
    # Paths
    if isinstance(id, float):
        id = f"{id:03.0f}"
    savename = figure_dir + f"/comparison_{id}"
    os.makedirs(savename, exist_ok=True)

    # Import data (lazy?)
    data_list = []
    for i in range(len(cases)):
        case = f"{cases[i]:02.0f}"
        print(f"Reading data for case {case}")
        data_list.append(xr.open_dataset(f"{output_dir}/data_repr_{case}.nc"))
    
    # Make figures
    vis_alongshore(data_list, title, cases, savename)
    vis_crossshore(data_list, title, cases, savename)

    return


### Executing Comparison
for i in range(num_comparisons):
    make_comparison(cases[i], titles[i])
