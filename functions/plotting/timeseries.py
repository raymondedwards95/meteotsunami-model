""" Functions for visualising time-series from model output

Main functions:
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib import gridspec

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
# fmt: on


class plot_timeseries():
    """ Methods to create a visualisation of time-series """
    number = 0

    def __init__(
        self,
        title: str = None,
    ):
        """ Create and setup a figure for time-series

        Options:
            title:  alternative title
        """
        plot_timeseries.number += 1

        self.figsize = FIGSIZE_NORMAL
        self.closed = False
        self.lines = np.array([0, 0, 0, 0])

        self.fig, self.axes = plt.subplots(4, 1)
        self.title = title

        print(
            f"\n# Initiated figure for time-series"
        )

    def _check_if_closed(self):
        """ Raises an error if the figure is supposed to be closed """
        if self.closed:
            raise TypeError(
                f"Figure is cleared and closed: it can not be edited")

    def _setup_figure(self):
        """ Figure setup """
        # Checks
        self._check_if_closed()
        if np.count_nonzero(self.lines) > 2:
            self.figsize = FIGSIZE_HIGH

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_tight_layout(True)

        # Figure specific
        if self.title is None:
            self.title = f"Time Series"
        self.fig.suptitle(self.title, va="top", ha="left", x=0.01)

    def _setup_plot(self):
        """ Plot setup """
        # Checks
        self._check_if_closed()

        # General
        for ax in self.axes:
            ax.axhline(color="black", linewidth=1)
            ax.legend(loc="upper right")
            ax.grid()

            max_lim = np.max(np.abs(ax.get_ylim()))
            ax.set_ylim(-1. * max_lim, max_lim)
            ax.set_xlim(0, None)
            ax.set_xlabel("Time [hours]")

        # wl
        ax = self.axes[0].set_ylabel("$wl$ [m]")

        # u
        self.axes[1].set_ylabel("$u$ [m/s]")

        # v
        self.axes[2].set_ylabel("$v$ [m/s]")

        # p
        self.axes[3].set_ylabel("$p$ [Pa]")
        self.axes[3].set_ylim(0, None)

    def _match_variable(self, variable: str) -> int:
        match variable.lower().strip():
            case "wl":
                return 0
            case "u":
                return 1
            case "v":
                return 2
            case "p":
                return 3
            case _:
                raise ValueError(
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

    def add_plot(
        self,
        dataset: xr.Dataset,
        variable: str,
        x: Numeric,
        y: Numeric,
        label: str = None,
    ):
        """ Adds data to a plot

        Input:
            dataset:    dataset containing gridded model output
            variable:   name of variable to plot
            x:          x-coordinate
            y:          y-coordinate

        Options:
            label:      label for plot
        """
        # Checks
        self._check_if_closed()

        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = "'unnamed'"

        ax_idx = self._match_variable(variable)
        self.lines[ax_idx] += 1

        if label is None:
            label = f"{dataset_name}"

        # Extract data
        time = dataset["t"] \
            .values \
            .astype("datetime64[s]") \
            .astype(float) / 3600.
        data = dataset[variable] \
            .interp(x=x, y=y)

        # Plot
        self.axes[ax_idx].plot(time, data, label=label)
        self.axes[ax_idx].fill_between(time, data, alpha=0.1)

        # Log
        print(
            f"# Added {variable} data from {dataset_name} for {x=} and {y=}, labeled by {label=}"
        )
        return self

    def save(
        self,
        saveloc: str,
        close: bool = True,
    ):
        """ Saves the figure

        Input:
            saveloc:    location where the figure should be saved

        Options:
            close:      close figure
        """
        # Checks
        self._check_if_closed()
        savename = f"{saveloc}/timeseries_{plot_timeseries.number:02.0f}" \
            .replace("//", "/")

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

        # update x-ticks
        self.axes[0].get_shared_x_axes() \
            .join(self.axes[0], *self.axes[1:])
        for ax in self.axes[:-1]:
            ax.set_xlabel("")
            ax.set_xticklabels([])

        # Save
        self.fig.execute_constrained_layout()
        self.fig.draw_without_rendering()
        self.fig.savefig(
            savename,
            bbox_inches="tight",
            dpi=FIG_DPI,
            pil_kwargs=FIG_PIL_KWARGS,
        )

        # End
        print(f"# Saved time-series figure as {savename}")
        if close:
            self.close()
        return

    def close(self):
        """ Close the figure """
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure time-series is closed")
        return


if __name__ == "__main__":
    # Additional imports
    import matplotlib.ticker as mticker
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    fmt = lambda x: f.format_data(x)

    # Define paths
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    figure_dir = f"{script_dir}/tests/figures/"
    os.makedirs(figure_dir, exist_ok=True)

    # Get data
    main_dir = os.path.dirname(script_dir)
    data_a = f"{main_dir}/reproduction-an-2012/output/data_repr_00.nc"
    data_a = xr.open_dataset(data_a, chunks="auto")

    # Make figures
    # fmt: off
    def test_timeseries():
        x = 1e4
        y = 1e6
        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "wl", x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .add_plot(data_a, "p",  x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "p",  x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .add_plot(data_a, "u",  x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "u",  x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .add_plot(data_a, "v",  x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "v",  x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=np.pi*y, label=f"$y={fmt(np.pi*y)}$") \
            .add_plot(data_a, "p",  x=x, y=np.pi*y, label=f"$y={fmt(np.pi*y)}$") \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "wl", x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .add_plot(data_a, "u",  x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "u",  x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .add_plot(data_a, "v",  x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "v",  x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .save(figure_dir)

        plot_timeseries() \
            .add_plot(data_a, "wl", x=x, y=1*y, label=f"$y={fmt(1*y)}$") \
            .add_plot(data_a, "wl", x=x, y=3*y, label=f"$y={fmt(3*y)}$") \
            .add_plot(data_a, "wl", x=x, y=5*y, label=f"$y={fmt(5*y)}$") \
            .save(figure_dir)

    # fmt: on
    test_timeseries()
