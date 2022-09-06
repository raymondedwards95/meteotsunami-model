""" Functions for visualising model outputs

Main functions:
"""

import os
import sys
import time
from typing import Union

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


class plot_spectrum_1d():
    """ Methods to create a visualisation of the 1d spectrum """

    def __init__(
        self,
        variable: str,
        demean: bool = True,
    ):
        """ Create and setup a figure for the 1d spectrum """
        self.figsize = FIGSIZE_NORMAL
        self.demean = demean

        self.fig, self.ax = plt.subplots(1, 1)

        self.variable: str
        self.variable_name: str
        match variable.lower().strip():
            case "wl":
                self.variable = "wl"
                self.variable_name = "Water level"
            case "u":
                self.variable = "u"
                self.variable_name = "x-component of water velocity"
            case "v":
                self.variable = "v"
                self.variable_name = "y-component of water velocity"
            case "p":
                self.variable = "p"
                self.variable_name = "Surface air pressure"
            case _:
                raise ValueError(
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'")

    def _setup_figure(
        self,
    ):
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_tight_layout(True)
        self.fig.suptitle(
            f"Power Spectrum - {self.variable_name}", va="top", ha="left", x=0.01)

    def _setup_plot(
        self,
    ):
        self.ax.axhline(color="black", linewidth=1)
        self.ax.set_ylim(0, None)
        self.ax.set_xlim(0, 2.2)
        self.ax.set_ylabel("Spectral Power [m$^2$ hr]")
        self.ax.legend(loc="upper right")
        return

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        y: Numeric,
        label: str = "",
    ):
        # Compute spectrum
        freqs, power = fa.spectral_analysis_1d(
            data=dataset,
            y=y,
            x=x,
            variable=self.variable,
            demean=self.demean,
        )

        # Scale time units from seconds to hour
        freqs *= 3600.
        power /= 3600.

        self.ax.plot(
            freqs,
            power,
            label=label
        )
        self.ax.fill_between(
            freqs,
            power,
            alpha=0.1,
        )
        return self

    def add_subplot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        y: Numeric,
    ):
        return self.add_plot(dataset, x, y)

    def save(
        self,
        saveloc: str,
    ):
        savename = f"{saveloc}_spectrum_1d"

        self._setup_plot()
        self.fig.savefig(
            savename,
            bbox_inches="tight",
            dpi=FIG_DPI,
            pil_kwargs=FIG_PIL_KWARGS,
        )
        return


if __name__ == "__main__":
    # Define paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    figure_dir = f"{script_dir}/tests/figures"
    figure_file = f"{figure_dir}/figure"
    os.makedirs(figure_dir, exist_ok=True)

    # Get data
    main_dir = os.path.dirname(script_dir)
    data_a = f"{main_dir}/reproduction-an-2012/output/data_repr_00.nc"
    data_b = f"{main_dir}/reproduction-an-2012/output/data_repr_01.nc"
    data_c = f"{main_dir}/reproduction-an-2012/output/data_repr_02.nc"

    data_a = xr.open_dataset(data_a, chunks="auto")
    data_b = xr.open_dataset(data_b, chunks="auto")
    data_c = xr.open_dataset(data_c, chunks="auto")

    # Make figures
    def test_spectrum_1d():
        x = 1e4
        y = 1e5
        var = "wl"
        plot_spectrum_1d() \
            .add_plot(data_a, var, x, y, label="a") \
            .add_plot(data_b, var, x, y, label="b") \
            .add_plot(data_c, var, x, y, label="c") \
            .save(figure_file)
    test_spectrum_1d()
