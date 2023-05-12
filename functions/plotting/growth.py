""" Functions for visualising wave evolution from model output

Main classes:
    plot_growth
"""

import os
import sys
from typing import Self

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
from functions.plotting.base import plot_base
# fmt: on


class plot_growth(plot_base):
    """Methods to create a visualisation of the evolution of wave height"""

    number = 0

    def __init__(
        self,
        title: str = None,
        x: Numeric = None,
    ) -> None:
        """Create and setup a figure for the evolution of wave height

        Options:
            `title`:    figure title
            `x`:        x-coordinate

        Methods:
            `add_plot`:     add data to the figure
            `save`:         write figure to disk as png and pgf
        """
        plot_growth.number += 1
        print("\n# Creating new figure")

        super().__init__()
        self.figsize = FIGSIZE_NORMAL

        self.fig, axes = plt.subplots(2, 1, sharex=True, squeeze=False)
        self.axes = axes.ravel()
        self.title = title
        self.figure_type = "Wave Height Evolution"
        self.figure_num = plot_growth.number
        self.line_num = -1

        self.time_max = 1

        self.x_ref: float
        if x is None:
            self.x_ref = 0
        else:
            self.x_ref = x

        self._check_if_closed()
        print(f"\n# Initiated figure '{self.figure_type} {self.figure_num}'")

    def _setup_figure(self) -> None:
        """Figure setup"""
        # Checks
        self._check_if_closed()

        # General
        super()._base_setup_figure()

    def _setup_plot(self) -> None:
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # General
        for ax in self.axes:
            ax.axhline(color="black", linewidth=1, alpha=0.5)
            ax.legend(loc="lower right")
            ax.grid()
            ax.set_xlim(0, self.time_max)
            ax.set_ylim(0, None)

        self.axes[1].set_xlabel("Time [\\si{\\hour}]")
        self.axes[0].set_ylabel("Water Level \n[\\si{\\meter}]")
        self.axes[1].set_ylabel("Surface Air Pressure \n[\\si{\\pascal}]")

        self.fig.align_labels()

    def add_plot(
        self,
        dataset: xr.Dataset,
        label: str = None,
    ) -> Self:
        """Adds data to a plot

        Input:
            `dataset`:  dataset containing gridded model output

        Options:
            `label`:    label for plot
        """
        # Checks
        self._check_if_closed()
        self.line_num += 1

        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = "-"

        if label is None:
            label = f"{dataset_name}"

        if np.isclose(self.x_ref, 0):
            self.x_ref = dataset["x"].min()

        # Plot data
        for j, var in enumerate(["wl", "p"]):
            data_time = dataset["t"].values.astype("datetime64[s]").astype(float) / 3600.0
            data_max = dataset[var].interp(x=self.x_ref).max(dim=["y"])
            data_min = -1.0 * dataset[var].interp(x=self.x_ref).min(dim=["y"])

            a = data_max > 0.99 * data_max.max()
            b = data_min > 0.99 * data_min.max()

            self.time_max = np.max(
                [
                    self.time_max,
                    data_time[a.argmax().values],
                    data_time[b.argmax().values],
                ]
            )

            self.axes[j].plot(
                data_time,
                data_max,
                color=f"C{self.line_num}",
                label=f"{label} - max",
            )
            self.axes[j].fill_between(
                data_time,
                data_max,
                color=f"C{self.line_num}",
                alpha=0.03,
            )

            self.axes[j].plot(
                data_time,
                data_min,
                color=f"C{self.line_num}",
                label=f"{label} - min",
                linestyle="--"
            )
            self.axes[j].fill_between(
                data_time,
                data_min,
                color=f"C{self.line_num}",
                alpha=0.03,
            )

        # Log
        print(
            f"# Added data from {dataset_name} at {self.x_ref}, labeled by {label=}"
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
        savename = f"{saveloc}/wavegrowth_{plot_growth.number:02.0f}".replace(
            "//", "/"
        )

        # Setup
        self._setup_figure()
        self._setup_plot()

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
    script_dir = f"{PATH_PLOTTING}"
    figure_dir = f"{PATH_TEST}"

    # Get data
    data_a = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_00.nc"
    data_b = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_01.nc"

    data_a = xr.open_dataset(data_a, chunks="auto")
    data_b = xr.open_dataset(data_b, chunks="auto")

    # Make figures
    # fmt: off
    def test_growth():
        plot_growth() \
            .add_plot(dataset=data_a) \
            .save(figure_dir)

        plot_growth() \
            .add_plot(dataset=data_a) \
            .add_plot(dataset=data_b) \
            .save(figure_dir)

        plot_growth(x=1e4) \
            .add_plot(dataset=data_a) \
            .add_plot(dataset=data_b) \
            .save(figure_dir)

        plot_growth(x=1e5) \
            .add_plot(dataset=data_a) \
            .add_plot(dataset=data_b) \
            .save(figure_dir)
    # fmt: on

    test_growth()
