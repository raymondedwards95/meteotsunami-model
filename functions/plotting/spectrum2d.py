""" Functions for visualising 2d spectra from model output

Main classes:
    plot_spectrum_1d
"""

import os
import sys
from typing import Self

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
import functions.analysis as fa
import functions.utilities as fu
from functions import *
from functions.plotting.base import plot_base
# fmt: on


class plot_spectrum_2d(plot_base):
    """Methods to create a visualisation of the 2d spectrum"""

    number = 0

    def __init__(
        self,
        variable: str,
        demean: bool = True,
        label: str = "",
        k_max: Numeric = None,
        f_max: Numeric = None,
        scale: str = "Mm",
        title: str = None,
    ) -> None:
        """Create and setup a figure for the 2d spectrum

        Input:
            `variable`:         name of variable, i.e. "wl", "u", "v" or "p"

        Options:
            `demean`:           remove mean from data
            `label`:            label for plot
            `k_max`:            upper limit for wavenumber k
            `f_max`:            upper limit for frequency f
            `scale`:            scale of plots ('m', 'km' or 'Mm')
            `title`:            figure title

        Methods:
            `add_plot`:         add data to the figure
            `add_dispersion`:   add theoretical relationship between wavenumber and frequency
            `save`:             writes the figure to disk
        """
        plot_spectrum_2d.number += 1
        self.time_scale_factor = 3600.0

        super().__init__()
        self.figsize = FIGSIZE_NORMAL

        self.fig, self.ax = plt.subplots(1, 1)
        self.title = title
        self.figure_type = "2D Power Spectrum"
        self.figure_num = plot_spectrum_2d.number
        self._setup_cax()

        self.demean = demean
        self.k_max = k_max
        self.f_max = f_max
        self.x = 0.0

        self.data_label = label

        self.unit = None
        self.scale_factor = None
        self.set_scale(scale)

        self.variable: str
        self.variable_long: str
        self.variable_unit: str
        self._set_variable(variable)

        self._check_if_closed()
        print(
            f"\n# Initiated figure '{self.figure_type} {self.figure_num}'",
            f"with variable '{self.variable}' ({self.variable_long.lower()})",
        )

    def _setup_figure(self) -> None:
        """Figure setup"""
        # Checks
        self._check_if_closed()

        # General
        if self.title is None:
            self.title = (
                f"Power Spectrum - {self.variable_long} - $x = \SI{{{self.x / 1000:0.1f}}}{{\kilo\meter}}$"
            )
        super()._base_setup_figure()

    def _setup_plot(self) -> None:
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # General
        self.ax.axhline(color="black", linewidth=1, alpha=0.5)
        self.ax.axvline(color="black", linewidth=1, alpha=0.5)
        self.ax.grid()
        # ax.ticklabel_format(scilimits=(-2, 2))

        # x-axis
        self.ax.set_xlabel(f"Wavenumber [1 / {self.unit}]")
        self.ax.set_xlim(0, fu.relative_ceil(0.8 * self.k_max_scaled, s=2))

        # y-axis
        self.ax.set_ylabel("Frequency [cycles per hour]")
        self.ax.set_ylim(0, fu.relative_ceil(0.8 * self.f_max_scaled, s=2))

    def _setup_cax(self) -> None:
        """Colorbar location setup

        Should be run after creating the figure and before adding data
        """
        # Checks
        self._check_if_closed()
        self._check_if_filled()

        # Make space for colorbar
        div = make_axes_locatable(self.ax)
        self.cax = div.append_axes(position="right", size="5%", pad="5%")

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        label: str = "",
        k_max: Numeric = None,
        f_max: Numeric = None,
        scale: str = None,
    ) -> Self:
        """Adds data to a plot

        Input:
            `dataset`:  dataset containing model output
            `x`:        x-coordinate

        Options:
            `label`:    label for plot
            `k_max`:    upper limit for wavenumber k
            `f_max`:    upper limit for frequency f
            `scale`:    scale of plots ('m', 'km' or 'Mm')
        """
        # Checks
        self._check_if_closed()
        self._check_if_filled()

        self.k_max = k_max
        self.f_max = f_max
        self.x = x

        if scale is not None:
            self.set_scale(scale)

        # Compute spectrum
        wavenumber, freqs, power = fa.spectral_analysis_2d(
            dataset,
            x=self.x,
            variable=self.variable,
            demean=self.demean,
        )
        # flip wavenumbers, since we only look at waves in a given direction
        wavenumber = -1.0 * wavenumber

        # Compute plot limits
        if self.k_max is None:
            self.k_max = 1.5 * np.max(wavenumber[np.argmax(power, axis=1)])

        if self.f_max is None:
            self.f_max = np.max(freqs[np.argmax(power, axis=0)])

        # Select cols and rows
        rows = freqs >= 0
        cols = wavenumber >= 0

        wavenumber = wavenumber[cols]
        freqs = freqs[rows]
        power = power[rows, :][:, cols]

        # Plot data
        self.im = self.ax.pcolormesh(
            wavenumber * self.scale_factor,
            freqs * self.time_scale_factor,
            power,
            shading="nearest",  # "gouraud",
            cmap=cmo.cm.matter,
            vmin=5.0 * np.min(power),
            rasterized=True,
        )

        # Colorbar
        self.cbar = self.fig.colorbar(self.im, cax=self.cax)
        self.cbar.set_label("Spectral power")  #TODO COMPUTE UNITS

        # Set plot limits
        self.k_max_scaled = np.min(
            [
                np.max(wavenumber) * self.scale_factor,
                fu.relative_ceil(self.k_max * self.scale_factor, s=2),
            ]
        )
        self.f_max_scaled = np.min(
            [
                np.max(freqs) * self.time_scale_factor,
                fu.relative_ceil(self.f_max * self.time_scale_factor, s=2),
            ]
        )

        # Log
        self.data_label = label
        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
            if self.data_label == "":
                self.data_label = dataset_name
        except KeyError:
            dataset_name = "-"
        print(f"# Added data from {dataset_name} for {x=}")
        return self

    def add_dispersion(
        self,
        n: Numeric,
        alpha: Numeric = 1.0 / 400.0,
        scale: str = None,
    ) -> Self:
        """Adds the dispersion relation to a plot

        Input:
            `n`:        maximum number of modes

        Options:
            `alpha`:    slope for computing dispersion relation
            `scale`:    scale of plots ('m', 'km' or 'Mm')
        """
        # Cheks
        self._check_if_closed()

        if scale is not None:
            self.set_scale(scale)

        # Get limits from figure
        xlims = self.ax.get_xlim()
        ylims = self.ax.get_ylim()

        # Compute and plot dispersion relation
        wavenumber = (
            np.linspace(np.min(xlims), self.k_max_scaled, 100) / self.scale_factor
        )

        for i in range(n):
            dispersion = fa.dispersion_relation(wavenumber, n=i, alpha=alpha)
            self.ax.plot(
                self.scale_factor * wavenumber,
                self.time_scale_factor * dispersion,
                linewidth=1,
                linestyle="--",
                color="black",
                rasterized=False,
            )
            self.ax.annotate(
                text=f"$n={i}$",
                xytext=(
                    0.8 * 0.99 * self.scale_factor * wavenumber[-1],
                    0.8 * 1.02 * self.time_scale_factor * dispersion[-1],
                ),
                xy=(
                    0.8 * self.scale_factor * wavenumber[-1],
                    0.8 * self.time_scale_factor * dispersion[-1],
                ),
                xycoords="data",
                textcoords="data",
                ha="right",
                va="top",
            )

        # Reset limits of figure
        self.ax.set_xlim(xlims)
        self.ax.set_ylim(ylims)

        # End
        print(f"# Added dispersion relation")
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
        savename = f"{saveloc}/spectrum_2d_{plot_spectrum_2d.number:02.0f}".replace(
            "//", "/"
        )

        # Setup
        self._setup_figure()
        self._setup_plot()

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
    def test_spectrum_2d():
        x = 1e4
        plot_spectrum_2d(variable="wl", demean=True) \
            .add_plot(data_a, x=x, scale="Mm") \
            .add_dispersion(n=3) \
            .save(figure_dir)

        plot_spectrum_2d(variable="wl", demean=True, scale="m") \
            .add_plot(data_b, x=x*2) \
            .add_dispersion(n=2) \
            .save(figure_dir)

        plot_spectrum_2d(variable="p", demean=False, scale="km") \
            .add_plot(data_a, x=x) \
            .save(figure_dir)

        plot_spectrum_2d(variable="u", demean=True) \
            .add_plot(data_a, x=x, scale="Mm") \
            .add_dispersion(n=2) \
            .save(figure_dir)

        plot_spectrum_2d(variable="v", demean=True) \
            .add_plot(data_a, x=x) \
            .add_dispersion(n=2) \
            .save(figure_dir)

        plot_spectrum_2d(variable="wl", demean=True) \
            .add_plot(data_a, x=x) \
            .add_dispersion(n=5) \
            .save(figure_dir)
    # fmt: on

    test_spectrum_2d()
