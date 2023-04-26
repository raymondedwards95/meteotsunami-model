""" Functions for visualising time-series from model output

Main classes:
    plot_timeseries
"""

import os
import sys
from typing import Self

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib import gridspec

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
from functions.plotting.base import plot_base
# fmt: on


class plot_timeseries(plot_base):
    """Methods to create a visualisation of time-series"""

    number = 0

    def __init__(
        self,
        title: str = None,
    ) -> None:
        """Create and setup a figure for time-series

        Options:
            `title`:    figure title

        Methods:
            `add_plot`: add data to the figure
            `save`:     write figure to disk as png and pgf
        """
        plot_timeseries.number += 1

        super().__init__()
        self.figsize = FIGSIZE_NORMAL
        self.lines = np.array([0, 0, 0, 0])

        self.fig, self.axes = plt.subplots(4, 1)
        self.title = title
        self.figure_type = "Time Series"
        self.figure_num = plot_timeseries.number

        self._check_if_closed()
        print(f"\n# Initiated figure '{self.figure_type} {self.figure_num}'")

    def _setup_figure(self) -> None:
        """Figure setup"""
        # Checks
        self._check_if_closed()
        if np.count_nonzero(self.lines) > 2:
            self.figsize = FIGSIZE_HIGH

        # General
        super()._base_setup_figure()

    def _setup_plot(self) -> None:
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # General
        for ax in self.axes:
            ax.axhline(color="black", linewidth=1, alpha=0.5)
            ax.legend(loc="upper right")
            ax.grid()
            # ax.ticklabel_format(scilimits=(-2, 2))

            max_lim = np.max(np.abs(ax.get_ylim()))
            ax.set_ylim(-1.0 * max_lim, max_lim)
            ax.set_xlim(0, None)
            ax.set_xlabel("Time [hours]")

        # wl
        ax = self.axes[0].set_ylabel("\\( wl \\) [\\si{\\meter}]")

        # u
        self.axes[1].set_ylabel("\\( u \\) [\\si{\\meter\\per\\second}]")

        # v
        self.axes[2].set_ylabel("\\( v \\) [\\si{\\meter\\per\\second}]")

        # p
        self.axes[3].set_ylabel("\\( p \\) [\\si{\\pascal}]")
        self.axes[3].set_ylim(0, None)

    def add_plot(
        self,
        dataset: xr.Dataset,
        variable: str,
        x: Numeric,
        y: Numeric,
        label: str = None,
    ) -> Self:
        """Adds data to a plot

        Input:
            `dataset`:  dataset containing gridded model output
            `variable`: name of variable to plot
            `x`:        x-coordinate
            `y`:        y-coordinate

        Options:
            `label`:    label for plot
        """
        # Checks
        self._check_if_closed()

        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = "-"

        ax_idx = self._match_variable(variable)
        self.lines[ax_idx] += 1

        if label is None:
            label = f"{dataset_name}; \\( x = \\SI{{{x / 1000.:0.1f}}}{{\\kilo\\meter}} \\); \\( y = \\SI{{{y/1000.:0.1f}}}{{\\kilo\\meter}} \\)"

        # Extract data
        time = dataset["t"].values.astype("datetime64[s]").astype(float) / 3600.0
        data = dataset[variable].interp(x=x, y=y)

        # Plot
        self.axes[ax_idx].plot(
            time,
            data,
            label=label,
            rasterized=False,
        )
        self.axes[ax_idx].fill_between(
            time,
            data,
            alpha=0.1,
            rasterized=False,
        )

        # Log
        print(
            f"# Added {variable} data from {dataset_name} for {x=} and {y=}, labeled by {label=}"
        )
        return self

    def save(
        self,
        saveloc: str,
        close: bool = True,
    ) -> None:
        """Saves the figure as png and pgf

        Input:
            `saveloc`:  location where the figure should be saved

        Options:
            `close`:    close figure
        """
        # Checks
        self._check_if_closed()
        savename = f"{saveloc}/timeseries_{plot_timeseries.number:02.0f}".replace(
            "//", "/"
        )

        # Setup
        self._setup_figure()
        self._setup_plot()

        # Remove unused axes/subplots
        for ax_idx in range(4):
            if self.lines[ax_idx] == 0:
                self.axes[ax_idx].remove()
        self.axes = self.axes[self.lines > 0]

        # Reposition axes/subplots
        gs = gridspec.GridSpec(nrows=len(self.axes), ncols=1, figure=self.fig)
        for ax_idx, ax in enumerate(self.axes):
            ax.set_subplotspec(gs[ax_idx])

        # update ticks
        for ax in self.axes:
            ax.label_outer()

        # Save
        self.fig.get_layout_engine().execute(self.fig)
        save_figure(self.fig, savename)

        # End
        print(f"# Saved figure '{self.figure_type} {self.figure_num}' as {savename}")
        if close:
            self.close()
        return


if __name__ == "__main__":
    # Additional imports
    import matplotlib.ticker as mticker

    f = mticker.ScalarFormatter(useOffset=False)

    def fmt(x):
        return f.format_data(x)

    # Define paths
    script_dir = f"{PATH_PLOTTING}"
    figure_dir = f"{PATH_TEST}"

    # Get data
    data_a = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_00.nc"
    data_b = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_01.nc"

    data_a = xr.open_dataset(data_a, chunks="auto")
    data_b = xr.open_dataset(data_b, chunks="auto")

    # Make figures
    # fmt: off
    def test_timeseries():
        x = 1e4
        y = 1e6
        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=1*y, label=f"\\( y = {fmt(1*y)} \\)") \
            .add_plot(data_a, "wl", x=x, y=2*y, label=f"\\( y = {fmt(2*y)} \\)") \
            .add_plot(data_a, "p",  x=x, y=1*y, label=f"\\( y = {fmt(1*y)} \\)") \
            .add_plot(data_a, "p",  x=x, y=2*y, label=f"\\( y = {fmt(2*y)} \\)") \
            .add_plot(data_a, "u",  x=x, y=1*y, label=f"\\( y = {fmt(1*y)} \\)") \
            .add_plot(data_a, "u",  x=x, y=2*y, label=f"\\( y = {fmt(2*y)} \\)") \
            .add_plot(data_a, "v",  x=x, y=1*y, label=f"\\( y = {fmt(1*y)} \\)") \
            .add_plot(data_a, "v",  x=x, y=2*y, label=f"\\( y = {fmt(2*y)} \\)") \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=np.pi*y) \
            .add_plot(data_a, "p",  x=x, y=np.pi*y) \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "u",  x=x, y=1*y, label=f"\\( y = {fmt(1*y)} \\)") \
            .add_plot(data_a, "u",  x=x, y=3*y, label=f"\\( y = {fmt(3*y)} \\)") \
            .add_plot(data_a, "u",  x=x, y=5*y, label=f"\\( y = {fmt(5*y)} \\)") \
            .add_plot(data_a, "v",  x=x, y=1*y, label=f"\\( y = {fmt(1*y)} \\)") \
            .add_plot(data_a, "v",  x=x, y=3*y, label=f"\\( y = {fmt(3*y)} \\)") \
            .add_plot(data_a, "v",  x=x, y=5*y, label=f"\\( y = {fmt(5*y)} \\)") \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=1*y) \
            .add_plot(data_a, "wl", x=x, y=3*y) \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=1*y) \
            .add_plot(data_a, "p", x=x, y=1*y) \
            .add_plot(data_a, "wl", x=x, y=3*y) \
            .add_plot(data_a, "p", x=x, y=3*y) \
            .add_plot(data_a, "wl", x=x, y=5*y) \
            .add_plot(data_a, "p", x=x, y=5*y) \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=np.pi*y) \
            .add_plot(data_a, "p",  x=x, y=np.pi*y) \
            .add_plot(data_b, "wl",  x=x, y=np.pi*y) \
            .add_plot(data_b, "p",  x=x, y=np.pi*y) \
            .save(figure_dir)
    # fmt: on

    test_timeseries()
