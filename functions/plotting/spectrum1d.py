""" Functions for visualising 1d spectra from model output

Main classes:
    plot_spectrum_1d
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
import functions.analysis as fa
from functions import *
from functions.plotting.base import plot_base
# fmt: on


class plot_spectrum_1d(plot_base):
    """Methods to create a visualisation of the 1d spectrum"""

    number = 0

    def __init__(
        self,
        variable: str,
        demean: bool = True,
        f_max: Numeric = None,
        title: str = None,
    ) -> None:
        """Create and setup a figure for the 1d spectrum

        Input:
            `variable`: name of variable, i.e. "wl", "u", "v" or "p"

        Options:
            `demean`:   remove mean from data
            `f_max`:    upper limit for frequency f (units: cycles per hour)
            `title`:    figure title

        Methods:
            `add_plot`: add data to the figure
            `save`:     writes the figure to disk
        """
        plot_spectrum_1d.number += 1

        super().__init__()
        self.figsize = FIGSIZE_NORMAL

        self.fig, self.ax = plt.subplots(1, 1)
        self.title = title
        self.figure_type = "1D Power Spectrum"
        self.figure_num = plot_spectrum_1d.number

        self.demean = demean
        self.f_max = f_max

        self.variable: str
        self.variable_long: str
        self.variable_unit: str
        self._set_variable(variable)

        self._check_if_closed()
        print(
            f"\n# Initiated figure '{self.figure_type} {self.figure_num}' with variable '{self.variable}' ({self.variable_long.lower()})"
        )

    def _setup_figure(self) -> None:
        """Figure setup"""
        # Checks
        self._check_if_closed()

        # General
        if self.title is None:
            self.title = f"Power Spectrum - {self.variable_long}"
        super()._base_setup_figure()

    def _setup_plot(self) -> None:
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # General
        self.ax.legend(loc="upper right")
        self.ax.grid()

        # x-axis
        self.ax.set_xlabel("Frequency [cycles per hour]")
        self.ax.set_xlim(0, self.f_max)
        # self.ax.xaxis.set_ticks(np.arange(0, self.f_max, 0.1))
        self.ax.xaxis.set_ticks(np.arange(0, self.f_max, 0.01), minor=True)

        # y-axis
        self.ax.set_ylabel(f"Spectral Power [{self.variable_unit}$^2$ hr]")
        self.ax.set_ylim(0, None)
        self.ax.ticklabel_format(
            axis="y",
            style="sci",
            scilimits=(0, 0),
            useMathText=True,
        )

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        y: Numeric,
        label: str = "",
    ) -> Self:
        """Adds data to a plot

        Input:
            `dataset`:  dataset containing model output
            `x`:        x-coordinate
            `y`:        y-coordinate

        Options:
            `label`:    label for plot
        """
        # Checks
        self._check_if_closed()

        # Compute spectrum
        freqs, power = fa.spectral_analysis_1d(
            data=dataset,
            y=y,
            x=x,
            variable=self.variable,
            demean=self.demean,
        )

        # Scale time units from seconds to hours
        freqs *= 3600.0
        power /= 3600.0

        # Set limits
        if self.f_max is None:
            self.f_max = 0.0

        new_f_max = freqs[  # take freq
            np.size(  # 'pick' last index
                np.trim_zeros(  # help finding last index
                    # find where power is significant
                    (power > np.var(power)).astype(int)
                )
            )
        ]

        if new_f_max > self.f_max:
            self.f_max = np.round(1.4 * new_f_max, 1)

        # Plot
        self.ax.plot(
            freqs,
            power,
            label=label,
            rasterized=False,
        )
        self.ax.fill_between(
            freqs,
            power,
            alpha=0.1,
            rasterized=False,
        )

        # Log
        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = "-"
        print(
            f"# Added data from {dataset_name} for {x=} and {y=}, labeled by {label=}"
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
        savename = f"{saveloc}/spectrum_1d_{plot_spectrum_1d.number:02.0f}".replace(
            "//", "/"
        )

        # Setup
        self._setup_figure()
        self._setup_plot()

        # Save
        self.fig.get_layout_engine().execute(self.fig)
        save_figure(self.fig, savename)

        # End
        print(
            f"# Saved figure '{self.figure_type} {self.figure_num}' figure as {savename}"
        )
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
    def test_spectrum_1d():
        x = 1e4
        y = 1e6
        plot_spectrum_1d(variable="wl", demean=True) \
            .add_plot(data_a, x, y, label="a") \
            .add_plot(data_b, x, y, label="b") \
            .save(figure_dir)

        plot_spectrum_1d(variable="wl", demean=False) \
            .add_plot(data_a, x, y, label=f"y={y}") \
            .add_plot(data_a, x, y*2, label=f"y={y*2}") \
            .add_plot(data_a, x, y*3, label=f"y={y*3}") \
            .save(figure_dir)

        plot_spectrum_1d(variable="wl", demean=False) \
            .add_plot(data_b, x, y, label=f"y={y}") \
            .add_plot(data_b, x, y*3, label=f"y={y*3}") \
            .add_plot(data_b, x, y*5, label=f"y={y*5}") \
            .add_plot(data_b, x, y*7, label=f"y={y*7}") \
            .save(figure_dir)
    # fmt: on

    test_spectrum_1d()
