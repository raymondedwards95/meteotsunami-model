""" Functions for visualising spectra from model output

Main classes:
    plot_spectrum_1d
    plot_spectrum_2d
"""

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
import functions.analysis as fa
import functions.utilities as fu
# fmt: on


class plot_spectrum_1d:
    """Methods to create a visualisation of the 1d spectrum"""

    number = 0

    def __init__(
        self,
        variable: str,
        demean: bool = True,
        f_max: Numeric = None,
    ):
        """Create and setup a figure for the 1d spectrum

        Input:
            `variable`: name of variable, i.e. "wl", "u", "v" or "p"

        Options:
            `demean`:   remove mean from data
            `f_max`:    upper limit for frequency f (units: cycles per hour)

        Methods:
            `add_plot`: add data to the figure
            `save`:     writes the figure to disk
        """
        plot_spectrum_1d.number += 1

        self.figsize = FIGSIZE_NORMAL
        self.closed = False
        self._check_if_closed()

        self.fig, self.ax = plt.subplots(1, 1)

        self.demean = demean
        self.f_max = f_max

        self.variable: str
        self.variable_long: str
        self.variable_unit: str
        self._set_variable(variable)

        print(
            f"\n# Initiated figure for spectrum_1d with variable '{self.variable}' ({self.variable_long.lower()})"
        )

    def _check_if_closed(self):
        """Raises an error if the figure is supposed to be closed"""
        if self.closed:
            raise TypeError(f"Figure is cleared and closed: it can not be edited")

    def _setup_figure(self):
        """Figure setup"""
        # Checks
        self._check_if_closed()

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_layout_engine("compressed")

        # Figure specific
        self.fig.suptitle(
            f"Power Spectrum - {self.variable_long}", va="top", ha="left", x=0.01
        )

    def _setup_plot(self):
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

    def _set_variable(self, variable: str):
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
                raise ValueError(f"{variable=} should be 'wl', 'u', 'v' or 'p'")

        return self

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        y: Numeric,
        label: str = "",
    ):
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
    ):
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
        print(f"# Saved 1d-spectrum figure as {savename}")
        if close:
            self.close()
        return

    def close(self):
        """Close the figure"""
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure 1d-spectrum is closed")
        return


class plot_spectrum_2d:
    """Methods to create a visualisation of the 2d spectrum"""

    number = 0

    def __init__(
        self,
        variable: str,
        demean: bool = True,
        label: str = "",
        k_max: Numeric = None,
        f_max: Numeric = None,
        scale="Mm",
    ):
        """Create and setup a figure for the 2d spectrum

        Input:
            `variable`:         name of variable, i.e. "wl", "u", "v" or "p"

        Options:
            `demean`:           remove mean from data
            `label`:            label for plot
            `k_max`:            upper limit for wavenumber k
            `f_max`:            upper limit for frequency f
            `scale`:            scale of plots ('m', 'km' or 'Mm')

        Methods:
            `add_plot`:         add data to the figure
            `add_dispersion`:   add theoretical relationship between wavenumber and frequency
            `save`:             writes the figure to disk
        """
        plot_spectrum_2d.number += 1
        self.y_scale = 3600.0

        self.figsize = FIGSIZE_NORMAL
        self.closed = False
        self.filled = False
        self._check_if_closed()

        self.fig, self.ax = plt.subplots(1, 1)
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

        print(
            f"\n# Initiated figure for spectrum_2d with variable '{self.variable}' ({self.variable_long.lower()})"
        )

    def _check_if_closed(self):
        """Raises an error if the figure is supposed to be closed"""
        if self.closed:
            raise TypeError(f"Figure is cleared and closed: it can not be edited")

    def _check_if_filled(self):
        """Raises an error if the figure is filled with data already"""
        if self.filled:
            raise TypeError(f"Figure contains data already: it cannot be edited")

    def _setup_figure(self):
        """Figure setup"""
        # Checks
        self._check_if_closed()

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_layout_engine("compressed")

        # Figure specific
        self.fig.suptitle(
            f"Power Spectrum - {self.variable_long} - $x = {self.x / 1000:0.1f}$ km",
            va="top",
            ha="left",
            x=0.01,
        )

    def _setup_plot(self):
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # General
        self.ax.axhline(color="black", linewidth=1, alpha=0.5)
        self.ax.axvline(color="black", linewidth=1, alpha=0.5)
        self.ax.grid()
        # ax.ticklabel_format(scilimits=(-2, 2), useMathText=True)

        # x-axis
        self.ax.set_xlabel(f"Wavenumber [1 / {self.unit}]")
        self.ax.set_xlim(0, fu.relative_ceil(self.k_max * self.scale_factor))

        # y-axis
        self.ax.set_ylabel("Frequency [cycles / hour]")
        self.ax.set_ylim(0, fu.relative_ceil(self.f_max * self.y_scale))

    def _setup_cax(self):
        """Colorbar location setup

        Should be run after creating the figure and before adding data
        """
        # Checks
        self._check_if_closed()
        self._check_if_filled()

        # Make space for colorbar
        div = make_axes_locatable(self.ax)
        self.cax = div.append_axes(position="right", size="5%", pad="5%")

    def _set_variable(self, variable: str):
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
                raise ValueError(f"{variable=} should be 'wl', 'u', 'v' or 'p'")

        return self

    def set_scale(self, scale: str):
        """Set the scale of all subplots

        Input:
            `scale`:    scale of plots ('m', 'km' or 'Mm')
        """
        self.unit: str
        self.scale_factor: float

        match scale:
            case "m":
                self.unit = "m"
                self.scale_factor = 1e0
            case "km":
                self.unit = "km"
                self.scale_factor = 1e3
            case "Mm":
                self.unit = "Mm"
                self.scale_factor = 1e6
            case _:
                raise ValueError(
                    f"Scale should be either 'm', 'km', or 'Mm', instead of '{scale}'"
                )

        return self

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        label: str = "",
        k_max: Numeric = None,
        f_max: Numeric = None,
        scale: str = None,
    ):
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
            freqs * self.y_scale,
            power,
            shading="nearest",  # "gouraud",
            cmap=cmo.cm.matter,
            vmin=5.0 * np.min(power),
            rasterized=True,
        )

        # Colorbar
        self.cbar = self.fig.colorbar(self.im, cax=self.cax)
        self.cbar.set_label("Spectral power")

        # Set plot limits
        self.ax.set_xlim(0, fu.relative_ceil(self.k_max * self.scale_factor))
        self.ax.set_ylim(0, fu.relative_ceil(self.f_max * self.y_scale))

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

    def add_dispersion(self, n: Numeric, alpha: Numeric = 1.0 / 400.0,
        scale: str = None,):
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
        wavenumber = np.linspace(np.min(xlims), np.max(xlims), 100) / self.scale_factor

        for i in range(n):
            dispersion = fa.dispersion_relation(wavenumber, n=i, alpha=alpha)
            self.ax.plot(
                self.scale_factor * wavenumber,
                self.y_scale * dispersion,
                linewidth=1,
                linestyle="--",
                color="black",
                rasterized=False,
            )
            self.ax.annotate(
                text=f"$n={i}$",
                xytext=(0.99 * self.scale_factor * wavenumber[-1], 1.02 * self.y_scale * dispersion[-1]),
                xy=(self.scale_factor * wavenumber[-1], self.y_scale * dispersion[-1]),
                xycoords="data",
                textcoords="data",
                ha="right",
                va="bottom",
            )

        # Reset limits of figure
        self.ax.set_xlim(xlims)
        self.ax.set_ylim(ylims)

        # End
        print(f"# Added dispersion relations to plot")
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
        print(f"# Saved 2d-spectrum figure as {savename}")
        if close:
            self.close()
        return

    def close(self):
        """Close the figure"""
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure 2d-spectrum is closed")
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

    def test_spectrum_2d():
        x = 1e4
        plot_spectrum_2d(variable="wl", demean=True) \
            .add_plot(data_a, x=x, scale="Mm") \
            .add_dispersion(n=3) \
            .save(figure_dir)

        plot_spectrum_2d(variable="wl", demean=True) \
            .add_plot(data_b, x=x*2) \
            .add_dispersion(n=2) \
            .save(figure_dir)

        plot_spectrum_2d(variable="p", demean=False, scale="km") \
            .add_plot(data_a, x=x) \
            .save(figure_dir)
    # fmt: on

    test_spectrum_1d()
    test_spectrum_2d()
