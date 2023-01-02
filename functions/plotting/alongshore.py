""" Functions for visualising along-shore wave profiles from model output

Main classes:
    plot_alongshore
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
from functions.plotting.base import plot_base
import functions.utilities as fu
# fmt: on


class plot_alongshore(plot_base):
    """Methods to create a visualisation of the along-shore profiles"""

    number = 0

    def __init__(
        self,
        variable: str,
        y_min: Numeric = None,
        y_max: Numeric = None,
        scale="Mm",
        title: str = None,
    ):
        """Create and setup a figure for along-shore profiles

        Input:
            `variable`:     name of variable to plot

        Options:
            `y_min`:        lower limit for y (in meters)
            `y_max`:        upper limit for y (in meters)
            `scale`:        scale of plots ('m', 'km' or 'Mm')
            `title`:    figure title

        Methods:
            `add_subplot`:  add a new subplot with data
            `add_plot`:     add data to the figure
        """
        plot_alongshore.number += 1

        super().__init__()
        self.figsize = FIGSIZE_NORMAL

        if y_min is None:
            self.y_min = 0
            self.y_min_fixed = False
        else:
            self.y_min = y_min
            self.y_min_fixed = True

        if y_max is None:
            self.y_max = 0
            self.y_max_fixed = False
        else:
            self.y_max = y_max
            self.y_max_fixed = True

        self.tx_done = set()

        self.fig = plt.figure()
        self.axes = []
        self.title = title
        self.figure_type = "Along Shore"
        self.figure_num = plot_alongshore.number

        self.unit = None
        self.scale_factor = None
        self.set_scale(scale)

        self.variable: str
        self.variable_long: str
        self.variable_unit: str
        self._set_variable(variable)

        self._check_if_closed()
        print(f"\n# Initiated figure '{self.figure_type} {self.figure_num}'")

    def _setup_figure(self):
        """Figure setup"""
        # Checks
        self._check_if_closed()
        if len(self.axes) > 2:
            self.figsize = FIGSIZE_HIGH
        if len(self.axes) > 4:
            self.figsize = FIGSIZE_LONG

        # General
        super()._base_setup_figure()

    def _setup_plot(self):
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # Compute y-limits
        max_lim = 0.0
        for ax in self.axes:
            ax_max_lim = np.max(np.abs(ax.get_ylim()))
            if ax_max_lim > max_lim:
                max_lim = ax_max_lim

        # General
        for ax in self.axes:
            ax.axhline(color="black", linewidth=1, alpha=0.5)
            ax.axvline(color="black", linewidth=1, alpha=0.5)
            ax.legend(loc="upper right")
            ax.grid()
            ax.set_xlim(
                fu.none_multiply(self.y_min, 1.0 / self.scale_factor),
                fu.none_multiply(self.y_max, 1.0 / self.scale_factor),
            )
            ax.set_ylim(-1.0 * max_lim, max_lim)
            # ax.ticklabel_format(scilimits=(-2, 2), useMathText=True)

        self.axes[-1].set_xlabel(f"$y$ [{self.unit}]")
        self.fig.supylabel(f"{self.variable_long} [{self.variable_unit}]")

    def add_subplot(
        self,
        dataset: xr.Dataset = None,
        x: Numeric = None,
        t: Numeric = None,
        label: str = None,
        y_min: Numeric = None,
        y_max: Numeric = None,
        scale: str = None,
    ):
        """Adds a subplot to the figure

        Can also add data to the new subplot using the following arguments:

        Input:
            `dataset`:  dataset containing gridded model output
            `x`:        x-coordinate
            `t`:        t-coordinate

        Options:
            `label`:    label for plot
            `y_min`:    lower limit for y (in kilometers)
            `y_max`:    upper limit for y (in kilometers)
            `scale`:    scale of plots ('m', 'km' or 'Mm')
        """
        # Checks
        self._check_if_closed()

        # Add subplot
        self.axes.append(self.fig.add_subplot())
        print(f"# Added subplot for '{self.figure_type} {self.figure_num}'")

        # Add data
        if (dataset is not None) and (x is not None) and (t is not None):
            return self.add_plot(
                dataset=dataset,
                x=x,
                t=t,
                label=label,
                y_min=y_min,
                y_max=y_max,
                scale=scale,
            )

        # End
        return self

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        t: Numeric,
        label: str = None,
        y_min: Numeric = None,
        y_max: Numeric = None,
        scale: str = None,
        skip_on_error: bool = False,
    ):
        """Adds data to a plot

        Input:
            `dataset`:          dataset containing gridded model output
            `x`:                x-coordinate (in meters)
            `t`:                t-coordinate

        Options:
            `label`:            label for plot
            `y_min`:            lower limit for y (in meters)
            `y_max`:            upper limit for y (in meters)
            `scale`:            scale of plots ('m', 'km' or 'Mm')
            `skip_on_error`:    skip plots when an error is catched
        """
        # Checks
        self._check_if_closed()

        if isinstance(x, (tuple, list, np.ndarray)) and np.size(x) > 1:
            print(f"# Given `{x=}` is a sequence")
            for x_single in x:
                self.add_plot(
                    dataset=dataset,
                    x=x_single,
                    t=t,
                    label=label,
                    y_min=y_min,
                    y_max=y_max,
                )
            return self

        if isinstance(t, (tuple, list, np.ndarray)) and np.size(t) > 1:
            print(f"# Given `{t=}` is a sequence")
            for t_single in t:
                self.add_plot(
                    dataset=dataset,
                    x=x,
                    t=t_single,
                    label=label,
                    y_min=y_min,
                    y_max=y_max,
                )
            return self

        if len(self.axes) == 0:
            self.add_subplot()

        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = "-"

        if label is None:
            label = f"{dataset_name}"

        if scale is not None:
            self.set_scale(scale)

        if x < dataset["x"].min():
            if skip_on_error:
                print(f"## {x=:0.1f} is outside the domain! Skipping...")
                return self
            print(
                f"## {x=:0.1f} is outside the domain!",
                f"Increasing x towards the domain.",
            )
            x = np.min(dataset["x"].values)
        if x > dataset["x"].max():
            if skip_on_error:
                print(f"## {x=:0.1f} is outside the domain! Skipping...")
                return self
            print(
                f"## {x=:0.1f} is outside the domain!",
                f"Decreasing x towards the domain.",
            )
            x = np.max(dataset["x"].values)

        if (x, t) in self.tx_done:
            print(f"## {x=:0.1f} and {t=:0.1f} is plotted before! Skipping...")
            return self
        self.tx_done.add((x, t))

        # Extract data
        data = dataset[self.variable].interp(x=x, t=fu.to_timestr(t))

        # Plot
        self.axes[-1].plot(
            dataset["y"] / self.scale_factor,
            data,
            label=label,
            rasterized=True,
        )
        self.axes[-1].fill_between(
            dataset["y"] / self.scale_factor,
            data,
            alpha=0.1,
            rasterized=True,
        )

        # Update plot-limits
        if self.y_min_fixed:
            pass
        else:
            if y_min is not None:  # new provided limits
                self.y_min_fixed = True
                self.y_min = y_min
            else:
                y_min_data = dataset["y"][np.abs(data) > 0.1 * data.std()][0]
                if self.y_min > (y_min_data):
                    self.y_min = y_min_data

        if self.y_max_fixed:
            pass
        else:
            if y_max is not None:  # new provided limits
                self.y_max_fixed = True
                self.y_max = y_max
            else:
                y_max_data = dataset["y"][np.abs(data) > 0.1 * data.std()][-1]
                if self.y_max < (y_max_data):
                    self.y_max = y_max_data

        # Log
        print(
            f"# Added {self.variable} data from {dataset_name}",
            f"for {x=} and {t=}, labeled by {label=}",
        )
        return self

    def save(
        self,
        saveloc: str,
        close: bool = True,
    ):
        """Saves the figure

        Input:
            `saveloc`:  location where the figure should be saved

        Options:
            `close`:    close figure
        """
        # Checks
        self._check_if_closed()
        savename = f"{saveloc}/along_{plot_alongshore.number:02.0f}".replace("//", "/")

        # Setup
        self._setup_figure()
        self._setup_plot()

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
    # Define paths
    figure_dir = f"{PATH_TEST}"

    # Get data
    data_a = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_00.nc"
    data_b = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_01.nc"

    data_a = xr.open_dataset(data_a, chunks="auto")
    data_b = xr.open_dataset(data_b, chunks="auto")

    # Make figures
    # fmt: off
    def test_alongshore():
        x = 1e4
        t = 36000
        plot_alongshore(variable="wl", scale="Mm") \
            .add_plot(data_a, x=x, t=t, label=f"a - t={t}") \
            .add_plot(data_a, x=x, t=2*t, label=f"a - t={2*t}") \
            .save(figure_dir)

        plot_alongshore(variable="wl") \
            .add_plot(data_a, x=[x, 2*x], t=t) \
            .save(figure_dir)

        plot_alongshore(variable="wl", y_max=5e6) \
            .add_subplot() \
            .add_plot(data_a, x=x, t=t, label=f"a - t={t}") \
            .add_subplot() \
            .add_plot(data_b, x=x, t=t, label=f"b - t={t}") \
            .save(figure_dir)

        plot_alongshore(variable="wl", y_min=0, scale="km") \
            .add_subplot(data_a, x=x, t=t, label=f"a - t={t}") \
            .add_subplot(data_a, x=x, t=2*t, label=f"a - t={2*t}") \
            .add_subplot(data_a, x=x, t=3*t, label=f"a - t={3*t}") \
            .add_subplot(data_a, x=x, t=4*t, label=f"a - t={4*t}") \
            .save(figure_dir)

        plot_alongshore(variable="wl", y_min=0) \
            .add_subplot(data_a, x=x, t=np.linspace(0, t, 5), scale="km") \
            .save(figure_dir)

        plot_alongshore(variable="wl", y_min=0) \
            .add_subplot(data_a, x=(x * np.arange(-5, 5, 2)), t=np.linspace(0, t, 5)) \
            .save(figure_dir)
    # fmt: on

    test_alongshore()
