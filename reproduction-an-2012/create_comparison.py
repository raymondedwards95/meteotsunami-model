"""Script for creating figures to compare with paper An et al. (2012) """

import argparse
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.utilities as fu
# fmt: on


def fig_1_contours_sse(dataset: xr.Dataset, saveloc: str) -> None:
    ta = time.perf_counter_ns()

    # Options
    savename = f"{saveloc}_01_sse_contours"
    t_list = np.array([4, 8, 12, 16]) * 1e4  # seconds
    y_min_list = np.array([0, 1, 2, 3]) * 1e3  # km
    y_max_list = np.array([4, 6, 8, 10]) * 1e3  # km

    # Figure
    print(f"Create figure '01_sse_contours'")
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=False)

    fig.supxlabel("$x$ [km]")
    fig.supylabel("$y$ [km]")

    # Subplots
    for i, ax in enumerate(np.ravel(axes)):
        t_single = t_list[i]
        print(f"# Add subplot {i}")

        ax.pcolormesh(
            dataset["x"] / 1e3,
            dataset["y"] / 1e3,
            dataset["wl"].interp(t=fu.to_timestr(t_single)),
            rasterized=True,
        )

        ax.set_xlim(0, 3e2)
        ax.set_ylim(y_min_list[i], y_max_list[i])

    # End
    tb = time.perf_counter_ns()
    save_figure(fig, savename)
    print(f"Saved figure as '{savename}'")
    print(f"Saved figure '01_sse_contours' in {(tb-ta)/1e6:0.0f} ms")
    return


def fig_2_along_sse(dataset: xr.Dataset, saveloc: str) -> None:
    ta = time.perf_counter_ns()

    # Options
    savename = f"{saveloc}_02_sse_along"
    t_list = np.array([4, 8, 12, 16]) * 1e4  # seconds
    y_min_list = np.array([0, 1, 2, 3]) * 1e3  # km
    y_max_list = np.array([4, 6, 8, 10]) * 1e3  # km
    x_single = 10e3  # m

    # Figure
    print(f"Create figure '02_sse_along'")
    fig, axes = plt.subplots(2, 2, sharex=False, sharey=True)

    fig.supxlabel("$x$ [km]")
    fig.supylabel("$y$ [km]")

    # Subplots
    for i, ax in enumerate(np.ravel(axes)):
        t_single = t_list[i]
        print(f"# Add subplot {i}")

        ax.plot(
            dataset["y"] / 1e3,
            dataset["wl"].interp(x=x_single, t=fu.to_timestr(t_single)),
        )

        ax.set_xlim(y_min_list[i], y_max_list[i])
        ax.set_ylim(-1, 1)

    # End
    tb = time.perf_counter_ns()
    save_figure(fig, savename)
    print(f"Saved figure as '{savename}'")
    print(f"Saved figure '02_sse_along' in {(tb-ta)/1e6:0.0f} ms")
    return


def fig_3_cross_sse(dataset: xr.Dataset, saveloc: str) -> None:
    ta = time.perf_counter_ns()

    # Options
    savename = f"{saveloc}_03_sse_cross"
    y_single = 7.56e6  # m; change this value to local max
    t_single = 1.6e5  # s

    # Figure
    print(f"Create figure '03_sse_cross'")
    fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)

    fig.supxlabel("$x$ [km]")
    fig.supylabel("$y$ [km]")

    # Subplots
    ax.plot(
        dataset["x"] / 1e3,
        dataset["wl"].interp(y=y_single, t=fu.to_timestr(t_single)),
    )

    ax.set_ylim(0, 0.8)
    ax.set_xlim(10, 600)

    # End
    tb = time.perf_counter_ns()
    save_figure(fig, savename)
    print(f"Saved figure as '{savename}'")
    print(f"Saved figure '03_sse_cross' in {(tb-ta)/1e6:0.0f} ms")
    return


if __name__ == "__main__":
    # Parameters
    t0 = time.perf_counter()

    parser = argparse.ArgumentParser(
        description="Create figures for a specific repr case to compare with paper An et al. (2012)",
    )
    parser.add_argument(
        "case",
        help="Case number",
        nargs="?",
        default="00",
        type=int,
    )
    args = parser.parse_args()

    case = f"{args.case:02}"

    # Paths
    current_dir = os.path.dirname(os.path.realpath(__file__))
    figure_dir = f"{current_dir}/compare/repr_{case}"
    data_dir = f"{current_dir}/output"

    data_file = f"{data_dir}/data_repr_{case}.nc"

    # Checks
    print(f"Making figures for repr {case=}")
    print(f"Save location: {figure_dir}")
    print(f"Data file:     {data_file}")

    if not os.path.exists(data_file):
        print(f"\nFile '{data_file}' is not found. Exiting.\n")
        sys.exit(1)

    # Get data
    data = xr.open_dataset(data_file, chunks={"y": -1, "x": -1, "t": "auto"})
    print(f"Grid size: {data.sizes}")
    print(f"Chunksize: {data.chunksizes}")

    # Make figures
    fig_1_contours_sse(dataset=data, saveloc=figure_dir)
    fig_2_along_sse(dataset=data, saveloc=figure_dir)
    fig_3_cross_sse(dataset=data, saveloc=figure_dir)
