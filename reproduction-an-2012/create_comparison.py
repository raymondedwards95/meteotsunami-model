""" Script to compare different simulations """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.analysis as fa
import functions.utilities as fu
import functions.visualisation as fv


### Cases
cases = [
    [0, 1, 2],
    [0, 3, 4, 5],
    [0, 10, 11, 12],
    [0, 15, 16, 17],
    [0, 20, 21],
    [0, 31, 32, 33],
    [0, 36],
    [36, 37, 38, 39]
]

titles = [
    "Maximum Computational Time Step",
    "Implicitness of Numerical Scheme",
    "Spatial Resolution - Pressure",
    "Temporal Resolution - Pressure",
    "Spatial Resolution - Model",
    "Pressure Distribution Location",
    "Bottom Slope",
    "Location Pressure and Slope Bottom"
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
    savename = savename.replace(".jpg", "") + "/along"
    _ylims = np.array([[0, 4], [1, 6], [2, 8], [3, 10]]) * 1e3
    _times = [4e4, 8e4, 12e4, 16e4]

    ## Figure
    fig, ax = plt.subplots(2, 2, sharex=False, sharey=True)
    fig.set_size_inches(8, 8)
    fig.set_dpi(150)
    fig.set_tight_layout(True)
    fig.suptitle(f"{title}\nAlong-shore Profile of Sea Surface Elevation")

    # Subplots
    for i in range(ax.size):
        _ax = ax[i//2, i % 2]  # select subplot
        _ax.grid()
        _ax.set_xlim(_ylims[i])
        _ax.set_ylim([-0.9, 0.9])
        _ax.set_title(f"$t = {_times[i]:0.0f}$s")
        _ax.axhline(color="black", linewidth=1)
        _ax.axvline(color="black", linewidth=1)
        if i//2: _ax.set_xlabel("$y$ [km]")
        if not i%2: _ax.set_ylabel("$SSE$ [m]")

        for j in range(len(data_list)):
            data = data_list[j]
            _ax.plot(
                data["y"]/1000.,
                data["wl"].interp(t=fu.to_timestr(_times[i]), x=10e3),
                color=f"C{j}",
                label=f"Case {cases[j]:02.0f}"
            )
        
        if i == 0:
            _ax.legend()
        
    fig.savefig(savename, bbox_inches="tight")
    print(f"Saved figure as '{savename}'")
    return

    
def vis_alongshore_diff(data_list, title, cases, savename):
    ## Parameters
    savename = savename.replace(".jpg", "") + "/along_diff"
    _ylims = np.array([[0, 4], [1, 6], [2, 8], [3, 10]]) * 1e3
    _times = [4e4, 8e4, 12e4, 16e4]

    ## Figure
    fig, ax = plt.subplots(2, 2, sharex=False, sharey=True)
    fig.set_size_inches(8, 8)
    fig.set_dpi(150)
    fig.set_tight_layout(True)
    fig.suptitle(f"{title}\nChange in Along-shore Profile")

    # Subplots
    for i in range(ax.size):
        _ax = ax[i//2, i % 2]  # select subplot
        _ax.grid()
        _ax.set_xlim(_ylims[i])
        _ax.set_ylim([-0.2, 0.2])
        _ax.set_title(f"$t = {_times[i]:0.0f}$s")
        _ax.axhline(color="black", linewidth=1)
        _ax.axvline(color="black", linewidth=1)
        if i//2: _ax.set_xlabel("$y$ [km]")
        if not i%2: _ax.set_ylabel("$\\Delta SSE$ [m]")

        reference = data_list[0]["wl"].interp(t=fu.to_timestr(_times[i]), x=10e3)

        for _j in range(len(data_list) - 1):
            j = _j + 1
            data = data_list[j]
            _ax.plot(
                data["y"]/1000.,
                data["wl"].interp(t=fu.to_timestr(_times[i]), x=10e3) - reference,
                color=f"C{j}",
                label=f"Case {cases[j]:02.0f}"
            )
        
        if i == 0:
            _ax.legend()
        
    fig.savefig(savename, bbox_inches="tight")
    print(f"Saved figure as '{savename}'")
    return


def vis_crossshore(data_list, title, cases, savename):
    ## Parameters
    savename = savename.replace(".jpg", "") + "/cross"
    _yslices = np.array([7.56]) * 1e6
    _tslice = 1.6e5

    ## Figure
    fig, ax = plt.subplots(_yslices.size, 1)
    fig.set_size_inches(8, 8)
    fig.set_dpi(150)
    fig.set_tight_layout(True)
    fig.suptitle(f"{title}\nCross-shore Profile of Sea Surface Elevation")
    
    if not isinstance(ax, np.ndarray):
        ax = np.array([ax])

    ## Subplots
    for i in range(ax.size):
        _ax = ax[i]
        _ax.axhline(color="black", linewidth=1)
        _ax.axvline(color="black", linewidth=1)
        _ax.grid()
        _ax.set_xlim([0, 600])
        _ax.set_ylim([0, 0.9])
        _ax.set_title(f"$y = {_yslices[i]/1000.}$km")
        _ax.set_ylabel("$SSE$ [m]")
        if i == ax.size: _ax.set_xlabel("x [km]")

        for j in range(len(data_list)):
            data = data_list[j]

            # Compute best fit parameters
            k0, y0 = fa.compute_decay_parameter(
                data, _yslices[i], _tslice
            )

            # Plot data
            plt.plot(
                data["x"] / 1000.,
                data["wl"].interp(t=fu.to_timestr(_tslice), y=_yslices[i]),
                color=f"C{j}",
                label=f"Case {cases[j]:02.0f}"
            )

            # Plot best fit
            plt.plot(
                data["x"] / 1000.,
                fa.exp_decay(data["x"], k0, y0),
                color=f"C{j}",
                linestyle="--",
                label=f"Best fit: $1/k_0 = {1./k0/1000.:0.1f}$km"
            )

        _ax.legend()

    fig.savefig(savename, bbox_inches="tight")
    print(f"Saved figure as '{savename}'")
    return


def vis_crossshore_diff(data_list, title, cases, savename):
    ## Parameters
    savename = savename.replace(".jpg", "") + "/cross_diff"
    _yslices = np.array([7.56]) * 1e6
    _tslice = 1.6e5

    ## Figure
    fig, ax = plt.subplots(_yslices.size, 1)
    fig.set_size_inches(8, 8)
    fig.set_dpi(150)
    fig.set_tight_layout(True)
    fig.suptitle(f"{title}\nChange in Cross-shore Profile")
    
    if not isinstance(ax, np.ndarray):
        ax = np.array([ax])

    ## Subplots
    for i in range(ax.size):
        _ax = ax[i]
        _ax.axhline(color="black", linewidth=1)
        _ax.axvline(color="black", linewidth=1)
        _ax.grid()
        _ax.set_xlim([0, 600])
        _ax.set_ylim([-0.2, 0.2])
        _ax.set_title(f"$y = {_yslices[i]/1000.}$km")
        _ax.set_ylabel("$\\Delta SSE$ [m]")
        if i == ax.size: _ax.set_xlabel("x [km]")

        reference = data_list[0]["wl"].interp(t=fu.to_timestr(_tslice), y=_yslices[i])

        for _j in range(len(data_list) - 1):
            j = _j + 1
            data = data_list[j]

            # Plot data
            plt.plot(
                data["x"] / 1000.,
                data["wl"].interp(t=fu.to_timestr(_tslice), y=_yslices[i]) - reference,
                color=f"C{j}",
                label=f"Case {cases[j]:02.0f}"
            )

        _ax.legend()

    fig.savefig(savename, bbox_inches="tight")
    print(f"Saved figure as '{savename}'")
    return


def make_comparison(cases, title, id="test"):
    # Paths
    if isinstance(id, float):
        id = f"{id:03.0f}"
    savename = figure_dir + f"/compare_{id}"
    os.makedirs(savename, exist_ok=True)

    # Import data (lazy?)
    data_list = []
    for i in range(len(cases)):
        case = f"{cases[i]:02.0f}"
        datafile = f"{output_dir}/data_repr_{case}.nc"
        print(f"Reading data for case {case} ({datafile})")
        data_list.append(xr.open_dataset(datafile))
    
    # Make figures
    vis_alongshore(data_list, title, cases, savename)
    vis_crossshore(data_list, title, cases, savename)
    vis_alongshore_diff(data_list, title, cases, savename)
    vis_crossshore_diff(data_list, title, cases, savename)

    return


### Executing Comparison
for i in range(num_comparisons):
    id = titles[i].lower().replace("-", "").replace("  ", " ").replace(" ", "_")
    id = f"{i:03.0f}_{id}"
    print(id)
    make_comparison(cases=cases[i], title=titles[i], id=id)

# make_comparison(cases=[0], title="test_1", id="test_1")
# make_comparison(cases=[0, 15], title="test_2", id="test_2")
# make_comparison(cases=[0, 15, 16], title="test_3", id="test_3")
# make_comparison(cases=[0, 15, 16, 17], title="test_4", id="test_4")
