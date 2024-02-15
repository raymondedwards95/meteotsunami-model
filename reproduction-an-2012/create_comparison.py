"""Script for creating figures to compare with paper An et al. (2012) """

import argparse
import os
import sys
import time

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.analysis as fa
import functions.utilities as fu
# fmt: on


# Helper functions
def closest_index_to_reference(array: npt.ArrayLike, value: Numeric) -> int:
    dist = np.abs(np.array(array) - value)
    idx = np.argmin(dist)
    return idx


# Figures
def fig_1_contours_sse(dataset: xr.Dataset, saveloc: str) -> None:
    ta = time.perf_counter_ns()

    # Options
    savename = f"{saveloc}_01_sse_contours"
    t_list = np.array([4, 8, 12, 16]) * 1e4  # seconds
    y_min_list = np.array([0, 1, 2, 3])  # Mm
    y_max_list = np.array([4, 6, 8, 10])  # Mm

    wl_absmax = 0.8  # np.nanmax(np.fabs(dataset["wl"]))

    subplot_labels = "abcd"

    # Figure
    print(f"Create figure '01_sse_contours'")
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=False)

    fig.set_size_inches(FIGSIZE_SMALL)
    fig.set_layout_engine("compressed")

    fig.supxlabel(
        "\\( x \\) [\\si{\\mega\\meter}]",
        x=0.8,
        y=0.0,
        va="top",
        ha="right",
        fontsize="medium",
    )
    fig.supylabel(
        "\\( y \\) [\\si{\\mega\\meter}]",
        x=0.0,
        y=0.975,
        va="top",
        ha="right",
        fontsize="medium",
    )

    # Subplots
    for i, ax in enumerate(np.ravel(axes)):
        t_single = t_list[i]
        print(f"# Add subplot {i}")

        im = ax.pcolormesh(
            dataset["x"] / 1e6,
            dataset["y"] / 1e6,
            dataset["wl"].interp(t=fu.to_timestr(t_single)),
            cmap=cmo.cm.balance,
            vmin=-1.0 * wl_absmax,
            vmax=wl_absmax,
            rasterized=True,
        )

        ax.annotate(
            f"({subplot_labels[i]})",
            (0.95, 0.95),
            xycoords="axes fraction",
            ha="right",
            va="top",
        )

        ax.set_xlim(0, 0.3)
        ax.set_ylim(y_min_list[i], y_max_list[i])

        ax.set_xticks(np.arange(0, 0.4, 0.1))
        ax.set_yticks(np.arange(y_min_list[i], y_max_list[i] + 1, 1))
        ax.ticklabel_format(scilimits=(-2, 2))

    # Colorbar
    cb = fig.colorbar(im, ax=axes, fraction=0.1, aspect=25, pad=0.01)
    cb.set_label("SSE [\\si{\\meter}]", loc="top", va="center")

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
    y_min_list = np.array([0, 1, 2, 3])  # Mm
    y_max_list = np.array([4, 6, 8, 10])  # Mm
    x_single = 10e3  # m

    subplot_labels = "abcd"

    # Figure
    print(f"Create figure '02_sse_along'")
    fig, axes = plt.subplots(2, 2, sharex=False, sharey=True)

    fig.set_size_inches(FIGSIZE_SMALL)
    fig.set_layout_engine("compressed")

    fig.supxlabel(
        "\\( y \\) [\\si{\\mega\\meter}]",
        x=0.975,
        y=0.0,
        va="top",
        ha="right",
        fontsize="medium",
    )
    fig.supylabel(
        "SSE [\\si{\\meter}]",
        x=0.0,
        y=0.975,
        va="top",
        ha="center",
        fontsize="medium",
    )

    # Subplots
    for i, ax in enumerate(np.ravel(axes)):
        t_single = t_list[i]
        print(f"# Add subplot {i}")
        wl = dataset["wl"].interp(x=x_single, t=fu.to_timestr(t_single))

        ax.plot(
            dataset["y"] / 1e6,
            wl,
            rasterized=True,
        )
        ax.fill_between(
            dataset["y"] / 1e6,
            wl,
            alpha=0.04,
            rasterized=True,
        )

        ax.axhline(
            color="black",
            linewidth=1,
            alpha=0.5,
        )

        ax.annotate(
            f"({subplot_labels[i]})",
            (0.95, 0.95),
            xycoords="axes fraction",
            ha="right",
            va="top",
        )

        ax.set_xlim(y_min_list[i], y_max_list[i])
        ax.set_ylim(-1, 1)

        ax.set_xticks(np.arange(y_min_list[i], y_max_list[i] + 1, 1))
        ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
        ax.ticklabel_format(scilimits=(-2, 2))

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

    # Find best y closest to given y_single
    y_idx = fu.find_local_maxima_y(
        dataset=dataset,
        t=t_single,
        x=10e3,
        variable="wl",
        minima=False,
        sort=True,
        number=5,
    )
    y_maxima = dataset["y"][y_idx]
    print(f"# Preferred value for y is {y_single:0.2e}")
    y_single = y_maxima[closest_index_to_reference(y_maxima, y_single)].values
    print(f"# Best value for y is {y_single:0.2e}")

    # Data
    wl = dataset["wl"].interp(y=y_single, t=fu.to_timestr(t_single))
    label_data = ""

    # Compute best fit
    k0, y0 = fa.compute_decay_parameter(
        dataset=dataset,
        y=y_single,
        t=t_single,
        variable="wl",
    )

    fit = fa.exp_decay(dataset["x"], k0, y0)
    label_fit = f"\\( A e^{{-k_0 x}} \\) with \\( 1 / k_0 = \\SI{{{1 / k0 / 1e3:0.0f}}}{{\\kilo\\meter}} \\)"

    # Figure
    print(f"Create figure '03_sse_cross'")
    fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)

    fig.set_size_inches(FIGSIZE_SMALL)
    fig.set_layout_engine("compressed")

    fig.supxlabel(
        "\\( x \\) [\\si{\\kilo\\meter}]",
        x=0.975,
        y=0.0,
        va="top",
        ha="right",
        fontsize="medium",
    )
    fig.supylabel(
        "SSE [\\si{\\meter}]",
        x=0.0,
        y=0.975,
        va="top",
        ha="right",
        fontsize="medium",
    )

    # Subplots
    ax.plot(
        dataset["x"] / 1e3,
        wl,
        linestyle="-",
        color="C0",
        label=label_data,
        rasterized=False,
    )
    ax.fill_between(  # TODO: maybe remove?
        dataset["x"] / 1e3,
        wl,
        color="C0",
        alpha=0.04,
        rasterized=False,
    )

    ax.plot(
        dataset["x"] / 1e3,
        fit,
        linestyle="--",
        color="C1",
        label=label_fit,
        rasterized=False,
    )
    ax.fill_between(  # TODO: maybe remove?
        dataset["x"] / 1e3,
        fit,
        color="C1",
        alpha=0.04,
        rasterized=False,
    )

    ax.set_ylim(0, 0.8)
    ax.set_xlim(10, 600)

    ax.set_xticks([10, 100, 200, 300, 400, 500, 600])
    ax.set_yticks(np.arange(0, 0.9, 0.1))

    ax.legend(fontsize="small")

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
