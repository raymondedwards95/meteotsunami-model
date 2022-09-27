""" Functions for visualising wave profiles from model output

Main classes:
    plot_alongshore
    plot_crossshore
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
import functions.utilities as fu
# fmt: on


class plot_alongshore():
    """ Methods to create a visualisation of the along-shore profiles """
    number = 0

    def __init__(
        self,
        variable: str,
    ):
        """ Create and setup a figure for along-shore profiles

        Input:
            `variable`:     name of variable to plot
        """
        plot_alongshore.number += 1

        self.figsize = FIGSIZE_NORMAL
        self.closed = False

        self.fig = plt.figure()
        self.axes = []
        self.title = None

        self.variable: str
        self.variable_long: str
        self.variable_unit: str
        match variable.lower().strip():
            case "wl":
                self.variable = "wl"
                self.variable_long = "Water level"
                self.variable_unit = "m"
            case "u":
                self.variable = "u"
                self.variable_long = "Cross shore water velocity"
                self.variable_unit = "m/s"
            case "v":
                self.variable = "v"
                self.variable_long = "Along shore water velocity"
                self.variable_unit = "m/s"
            case "p":
                self.variable = "p"
                self.variable_long = "Surface air pressure"
                self.variable_unit = "Pa"
            case _:
                raise ValueError(
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

        print(
            f"\n# Initiated figure for along-shore profiles"
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
        if len(self.axes) > 2:
            self.figsize = FIGSIZE_HIGH
        if len(self.axes) > 4:
            self.figsize = FIGSIZE_LONG

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_tight_layout(True)

        # Figure specific
        if self.title is None:
            self.title = f"Along-shore profile"
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

            ax.set_xlabel(f"$y$ [km]")

        self.axes[len(self.axes) // 2].set_ylabel(f"{self.variable_long} [{self.variable_unit}]")

    def add_subplot(
        self,
        dataset: xr.Dataset = None,
        x: Numeric = None,
        t: Numeric = None,
        label: str = None,
    ):
        """ Adds a subplot to the figure

        Can also add data to the new subplot using the following arguments:

        Input:
            `dataset`:  dataset containing gridded model output
            `variable`: name of variable to plot
            `x`:        x-coordinate
            `t`:        t-coordinate

        Options:
            `label`:    label for plot
        """
        # Checks
        self._check_if_closed()

        # Add subplot
        self.axes.append(self.fig.add_subplot())
        print(f"# Added subplot for along-shore profiles")

        # Add data
        if (dataset is not None) and (x is not None) and (t is not None):
            return self.add_plot(
                dataset=dataset,
                x=x,
                t=t,
                label=label,
            )

        # End
        return self

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        t: Numeric,
        label: str = None,
    ):
        """ Adds data to a plot

        Input:
            `dataset`:  dataset containing gridded model output
            `variable`: name of variable to plot
            `x`:        x-coordinate
            `t`:        t-coordinate

        Options:
            `label`:    label for plot
        """
        # Checks
        self._check_if_closed()

        if len(self.axes) == 0:
            self.add_subplot()

        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = "'unnamed'"

        if label is None:
            label = f"{dataset_name}"

        # Extract data
        y = dataset["y"] / 1000.
        data = dataset[self.variable] \
            .interp(x=x, t=fu.to_timestr(t))

        # Plot
        self.axes[-1].plot(y, data, label=label)
        self.axes[-1].fill_between(y, data, alpha=0.1)

        # Log
        print(
            f"# Added {self.variable} data from {dataset_name} for {x=} and {t=}, labeled by {label=}"
        )
        return self

    def save(
        self,
        saveloc: str,
        close: bool = True,
    ):
        """ Saves the figure

        Input:
            `saveloc`:  location where the figure should be saved

        Options:
            `close`:    close figure
        """
        # Checks
        self._check_if_closed()
        savename = f"{saveloc}/along_{plot_alongshore.number:02.0f}" \
            .replace("//", "/")

        # Setup
        self._setup_figure()
        self._setup_plot()

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
        self.fig.savefig(
            savename,
            bbox_inches="tight",
            dpi=FIG_DPI,
            pil_kwargs=FIG_PIL_KWARGS,
        )

        # End
        print(f"# Saved along-shore figure as {savename}")
        if close:
            self.close()
        return

    def close(self):
        """ Close the figure """
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure along-shore is closed")
        return


class plot_crossshore():
    """ Methods to create a visualisation of the cross-shore profiles """
    number = 0

    def __init__(
        self,
        variable: str,
    ):
        """ Create and setup a figure for cross-shore profiles

        Input:
            `variable`:     name of variable to plot
        """
        plot_crossshore.number += 1

        self.figsize = FIGSIZE_NORMAL
        self.closed = False

        self.fig = plt.figure()

        self.variable: str
        match variable.lower().strip():
            case "wl":
                self.variable = "wl"
            case "u":
                self.variable = "u"
            case "v":
                self.variable = "v"
            case "p":
                self.variable = "p"
            case _:
                raise ValueError(
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )


if __name__ == "__main__":
    # Define paths
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    figure_dir = f"{script_dir}/tests/figures/"
    os.makedirs(figure_dir, exist_ok=True)

    # Get data
    main_dir = os.path.dirname(script_dir)
    data_a = f"{main_dir}/reproduction-an-2012/output/data_repr_00.nc"
    data_b = f"{main_dir}/reproduction-an-2012/output/data_repr_01.nc"

    data_a = xr.open_dataset(data_a, chunks="auto")
    data_b = xr.open_dataset(data_b, chunks="auto")

    # Make figures
    # fmt: off
    def test_alongshore():
        x = 1e4
        t = 3600
        plot_alongshore(variable="wl") \
            .add_plot(data_a, x=x, t=t, label=f"a - t={t}") \
            .add_plot(data_a, x=x, t=2*t, label=f"a - t={2*t}") \
            .save(figure_dir)

        plot_alongshore(variable="wl") \
            .add_subplot() \
            .add_plot(data_a, x=x, t=t, label=f"a - t={t}") \
            .add_subplot() \
            .add_plot(data_b, x=x, t=t, label=f"b - t={t}") \
            .save(figure_dir)

        plot_alongshore(variable="wl") \
            .add_subplot(data_a, x=x, t=t, label=f"a - t={t}") \
            .add_subplot(data_a, x=x, t=2*t, label=f"a - t={2*t}") \
            .add_subplot(data_a, x=x, t=3*t, label=f"a - t={3*t}") \
            .add_subplot(data_a, x=x, t=4*t, label=f"a - t={4*t}") \
            .add_subplot(data_a, x=x, t=5*t, label=f"a - t={5*t}") \
            .save(figure_dir)

    # fmt: on
    test_alongshore()
