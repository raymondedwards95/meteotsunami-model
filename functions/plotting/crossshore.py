""" Functions for visualising cross-shore wave profiles from model output

Main classes:
    plot_crossshore
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
import functions.analysis as fa
import functions.utilities as fu
from functions import *
from functions.plotting.base import plot_base
# fmt: on


class plot_crossshore(plot_base):
    """Methods to create a visualisation of the cross-shore profiles"""

    number = 0

    def __init__(
        self,
        variable: str,
        x_max: Numeric = None,
        scale="Mm",
        title: str = None,
    ) -> None:
        """Create and setup a figure for cross-shore profiles

        Input:
            `variable`:     name of variable to plot

        Options:
            `x_max`:        upper limit for x (in meters)
            `scale`:        scale of plots ('m', 'km' or 'Mm')
            `title`:    figure title

        Methods:
            `plot_peaks`:   plot cross-shore profiles of peaks in different subplots
            `add_subplot`:  add a new subplot with data
            `add_plot`:     add data to the figure
        """
        plot_crossshore.number += 1

        super().__init__()
        self.figsize = FIGSIZE_NORMAL

        if x_max is None:
            self.x_max = 0
            self.x_max_fixed = False
        else:
            self.x_max = x_max
            self.x_max_fixed = True

        self.fig = plt.figure()
        self.axes = []
        self.title = title
        self.figure_type = "Cross Shore"
        self.figure_num = plot_crossshore.number

        self.unit = None
        self.scale_factor = None
        self.set_scale(scale)

        self.variable: str
        self.variable_long: str
        self.variable_unit: str
        self._set_variable(variable)

        self._check_if_closed()
        print(f"\n# Initiated figure '{self.figure_type} {self.figure_num}'")

    def _setup_figure(self) -> None:
        """Figure setup"""
        # Checks
        self._check_if_closed()
        if len(self.axes) > 2:
            self.figsize = FIGSIZE_HIGH
        if len(self.axes) > 4:
            self.figsize = FIGSIZE_LONG

        # General
        super()._base_setup_figure()

    def _setup_plot(self) -> None:
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # Compute y-limits
        max_lim = 0.0
        min_lim = 0.0
        for ax in self.axes:
            ax_max_lim = np.max(ax.get_ylim())
            ax_min_lim = np.min(ax.get_ylim())

            if ax_max_lim > max_lim:
                max_lim = ax_max_lim

            if ax_min_lim < min_lim:
                min_lim = ax_min_lim

        # General
        for ax in self.axes:
            ax.axhline(color="black", linewidth=1, alpha=0.5)
            ax.axvline(color="black", linewidth=1, alpha=0.5)
            ax.legend(loc="upper right")
            ax.grid()
            ax.set_xlim(0, fu.none_multiply(self.x_max, 1.0 / self.scale_factor))
            ax.set_ylim(min_lim, max_lim)
            ax.ticklabel_format(scilimits=(-3, 3))

        self.axes[-1].set_xlabel(f"$x$ [{self.unit}]")
        self.fig.supylabel(f"{self.variable_long} [{self.variable_unit}]")

    def plot_peaks(
        self,
        dataset: xr.Dataset,
        t: Numeric,
        fit_curve: bool = True,
        x_max: Numeric = None,
        sort: bool = False,
        number: Integer = None,
        scale: str = None,
    ) -> Self:
        """Finds local maxima for fixed t and plots the cross-shore profiles in separate subplots

        Input:
            `dataset`:      dataset containing gridded model output
            `t`:            t-coordinate

        Options:
            `fit_curve`:    apply a curve fit on the data
            `x_max`:        upper limit for x (in meters)
            `sort`:         sort plots by highest values in data
            `number`:       number of waves to plot
            `scale`:        scale of plots ('m', 'km' or 'Mm')
        """
        # Checks
        self._check_if_closed()

        # Find indices and y-coordinates
        x = np.min([dataset["x"][10], dataset["x"].max() / 20.0])
        y_idx = fu.find_local_maxima_y(
            dataset,
            t,
            x,
            variable=self.variable,
            minima=False,
            sort=sort,
            number=number,
        )

        if len(y_idx) < 1:
            print(f"# There are no local maxima to plot!")
            self.axes.append(self.fig.add_subplot())  # create empty plot
            return self

        # Make plots
        for idx in y_idx:
            y = dataset["y"][idx].values
            self.add_subplot(
                dataset=dataset,
                t=t,
                y=y,
                fit_curve=fit_curve,
                label=None,
                x_max=x_max,
                scale=scale,
            )

        # End
        print(f"# Added {len(y_idx)} plots")
        return self

    def add_subplot(
        self,
        dataset: xr.Dataset = None,
        t: Numeric = None,
        y: Numeric = None,
        fit_curve: bool = True,
        label: str = None,
        x_max: Numeric = None,
        scale: str = None,
    ) -> Self:
        """Adds a subplot to the figure

        Can also add data to the new subplot using the following arguments:

        Input:
            `dataset`:      dataset containing gridded model output
            `t`:            t-coordinate
            `y`:            y-coordinate

        Options:
            `fit_curve`:    apply a curve fit on the data
            `label`:        label for plot
            `x_max`:        upper limit for x (in meters)
            `scale`:        scale of plots ('m', 'km' or 'Mm')
        """
        # Checks
        self._check_if_closed()

        # Add subplot
        self.axes.append(self.fig.add_subplot())
        print(f"# Added subplot for '{self.figure_type} {self.figure_num}'")

        # Add data
        if (dataset is not None) and (y is not None) and (t is not None):
            return self.add_plot(
                dataset=dataset,
                t=t,
                y=y,
                fit_curve=fit_curve,
                label=label,
                x_max=x_max,
                scale=scale,
            )

        # End
        return self

    def add_plot(
        self,
        dataset: xr.Dataset,
        t: Numeric,
        y: Numeric,
        fit_curve: bool = True,
        label: str = None,
        x_max: Numeric = None,
        scale: str = None,
    ) -> Self:
        """Adds data to a plot

        Input:
            `dataset`:      dataset containing gridded model output
            `t`:            t-coordinate
            `y`:            y-coordinate

        Options:
            `fit_curve`:    apply a curve fit on the data
            `label`:        label for plot
            `x_max`:        upper limit for x (in meters)
            `scale`:        scale of plots ('m', 'km' or 'Mm')
        """
        # Checks
        self._check_if_closed()

        if isinstance(y, (tuple, list, np.ndarray)) and np.size(y) > 1:
            print(f"# Given `{y=}` is a sequence")
            for y_single in y:
                self.add_plot(
                    dataset=dataset,
                    t=t,
                    y=y_single,
                    fit_curve=fit_curve,
                    label=label,
                    x_max=x_max,
                    scale=scale,
                )
            return self

        if isinstance(t, (tuple, list, np.ndarray)) and np.size(t) > 1:
            print(f"# Given `{t=}` is a sequence")
            for t_single in t:
                self.add_plot(
                    dataset=dataset,
                    t=t_single,
                    y=y,
                    fit_curve=fit_curve,
                    label=label,
                    x_max=x_max,
                    scale=scale,
                )
            return self

        if len(self.axes) == 0:
            self.add_subplot()

        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = ""

        if label is None:
            label = f"$y = \\SI{{{y/1e3:0.0f}}}{{\\kilo\\meter}}$, $t = \\SI{{{t/3600.:0.1f}}}{{\\hour}}$"

        if scale is not None:
            self.set_scale(scale)

        # Extract data
        data = dataset[self.variable].interp(y=y, t=fu.to_timestr(t))

        # Plot
        self.axes[-1].plot(
            dataset["x"] / self.scale_factor,
            data,
            label=label,
            rasterized=False,
        )
        self.axes[-1].fill_between(
            dataset["x"] / self.scale_factor,
            data,
            alpha=0.1,
            rasterized=False,
        )

        # Fit curve
        if fit_curve:
            k0, y0 = fa.compute_decay_parameter(
                dataset=dataset, y=y, t=t, variable=self.variable
            )

            fit = fa.exp_decay(dataset["x"], k0, y0)
            label_fit = f"$A e^{{-k_0 x}}$ with $1/k_0 = \\SI{{{1/k0/1e3:0.0f}}}{{\\kilo\\meter}}$"

            self.axes[-1].plot(
                dataset["x"] / self.scale_factor,
                fit,
                "--",
                label=label_fit,
                rasterized=False,
            )
            # self.axes[-1].fill_between(dataset["x"] / self.scale_factor, fit, alpha=0.1, rasterized=False,)

        # Update plot-limits
        if self.x_max_fixed:
            pass
        else:
            if x_max is not None:  # new provided limits
                self.x_max_fixed = True
                self.x_max = x_max
            else:
                x_max_data = dataset["x"].max()
                if self.x_max < (x_max_data):
                    self.x_max = x_max_data

        # Log
        print(
            f"# Added {self.variable} data from {dataset_name} for {y=} and {t=}, labeled by {label=}"
        )
        return self

    def save(
        self,
        saveloc: str,
        close: bool = True,
    ) -> None:
        """Saves the figure

        Input:
            `saveloc`:  location where the figure should be saved

        Options:
            `close`:    close figure
        """
        # Checks
        self._check_if_closed()
        savename = f"{saveloc}/cross_{plot_crossshore.number:02.0f}".replace("//", "/")

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
    def test_crossshore():
        y = [52e5, 72e5]
        t = 3600*42
        plot_crossshore(variable="wl", scale="km") \
            .add_plot(data_a, t=t, y=y[0], label=f"a - y={y[0]}") \
            .add_plot(data_a, t=t, y=y[1], label=f"a - y={y[1]}") \
            .save(figure_dir)

        plot_crossshore(variable="wl", x_max=5e5) \
            .add_subplot(data_a, t=t, y=y, scale="km") \
            .save(figure_dir)

        plot_crossshore(variable="wl", x_max=3e5) \
            .plot_peaks(data_a, t) \
            .save(figure_dir)

        plot_crossshore(variable="wl", x_max=5e5) \
            .plot_peaks(data_a, t, number=3, scale="km") \
            .save(figure_dir)

        plot_crossshore(variable="wl", x_max=7e5) \
            .plot_peaks(data_a, t, number=7, scale="m") \
            .save(figure_dir)

        plot_crossshore(variable="wl", x_max=3.5e5) \
            .plot_peaks(data_a, t, number=3, sort=True) \
            .save(figure_dir)
    # fmt: on

    test_crossshore()
