""" Functions for visualising spectra from model output

Main functions:
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
# fmt: on


class plot_spectrum_1d():
    """ Methods to create a visualisation of the 1d spectrum """
    number = 0

    def __init__(
        self,
        variable: str,
        demean: bool = True,
    ):
        """ Create and setup a figure for the 1d spectrum

        Input:
            variable:   name of variable, i.e. "wl", "u", "v" or "p"

        Options:
            demean:     remove mean from data
        """
        plot_spectrum_1d.number += 1

        self.figsize = FIGSIZE_NORMAL
        self.demean = demean
        self.closed = False
        self._check_if_closed()

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
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

        print(
            f"\n# Initiated figure for spectrum_1d with variable '{self.variable}' ({self.variable_name.lower()})"
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

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_tight_layout(True)

        # Figure specific
        self.fig.suptitle(
            f"Power Spectrum - {self.variable_name}", va="top", ha="left", x=0.01
        )

    def _setup_plot(self):
        """ Plot setup """
        # Checks
        self._check_if_closed()

        # General
        self.ax.axhline(color="black", linewidth=1)
        self.ax.legend(loc="upper right")
        self.ax.grid()

        # x-axis
        self.ax.set_xlabel("Frequency [cycles per hour]")
        self.ax.set_xlim(0, 2.2)
        self.ax.xaxis.set_ticks(np.arange(0, 2.2, 0.2))
        self.ax.xaxis.set_ticks(np.arange(0, 2.2, 0.1), minor=True)

        # y-axis
        self.ax.set_ylabel("Spectral Power [m$^2$ hr]")
        self.ax.set_ylim(0, None)
        self.ax.ticklabel_format(
            axis="y", style="sci", scilimits=(0, 0), useMathText=True,
        )

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        y: Numeric,
        label: str = "",
    ):
        """ Adds data to a plot

        Input:
            dataset:    dataset containing model output
            x:          x-coordinate
            y:          y-coordinate

        Options:
            label:      label for plot
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
        freqs *= 3600.
        power /= 3600.

        # Plot
        self.ax.plot(
            freqs,
            power,
            label=label,
        )
        self.ax.fill_between(
            freqs,
            power,
            alpha=0.1,
        )

        # Log
        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
        except KeyError:
            dataset_name = "'unnamed'"
        print(
            f"# Added data from {dataset_name} for {x=} and {y=}, labeled by {label=}"
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
        savename = f"{saveloc}/spectrum_1d_{plot_spectrum_1d.number:02.0f}" \
            .replace("//", "/")

        # Setup
        self._setup_figure()
        self._setup_plot()

        # Save
        self.fig.savefig(
            savename,
            bbox_inches="tight",
            dpi=FIG_DPI,
            pil_kwargs=FIG_PIL_KWARGS,
        )

        # End
        print(f"# Saved 1d-spectrum figure as {savename}")
        if close:
            self.close()
        return

    def close(self):
        """ Close the figure """
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure 1d-spectrum is closed")
        return


class plot_spectrum_2d():
    """ Methods to create a visualisation of the 2d spectrum """
    number = 0

    def __init__(
        self,
        variable: str,
        demean: bool = True,
        label: str = "",
    ):
        """ Create and setup a figure for the 2d spectrum

        Input:
            variable:   name of variable, i.e. "wl", "u", "v" or "p"

        Options:
            demean:     remove mean from data
            label:      label for plot
        """
        plot_spectrum_2d.number += 1
        self.x_scale = 1e6
        self.y_scale = 3600.

        self.figsize = FIGSIZE_NORMAL
        self.demean = demean
        self.closed = False
        self.filled = False
        self._check_if_closed()

        self.fig, self.ax = plt.subplots(1, 1)
        self._setup_cax()

        self.x = 0.

        self.data_label = label

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
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

        print(
            f"\n# Initiated figure for spectrum_2d with variable '{self.variable}' ({self.variable_name.lower()})"
        )

    def _check_if_closed(self):
        """ Raises an error if the figure is supposed to be closed """
        if self.closed:
            raise TypeError(
                f"Figure is cleared and closed: it can not be edited")

    def _check_if_filled(self):
        """ Raises an error if the figure is filled with data already """
        if self.filled:
            raise TypeError(
                f"Figure contains data already: it cannot be edited")

    def _setup_figure(self):
        """ Figure setup """
        # Checks
        self._check_if_closed()

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_tight_layout(True)

        # Figure specific
        self.fig.suptitle(
            f"Power Spectrum - {self.variable_name} - {self.data_label} - $x = {self.x / 1000:0.1f}$ km",
            va="top",
            ha="left",
            x=0.01
        )

    def _setup_plot(self):
        """ Plot setup """
        # Checks
        self._check_if_closed()

        # General
        self.ax.axhline(color="black", linewidth=1)
        self.ax.axvline(color="black", linewidth=1)
        self.ax.grid()

        # x-axis
        self.ax.set_xlabel("Wavenumber [1 / (1000 km)]")

        # y-axis
        self.ax.set_ylabel("Frequency [cycles / hour]")
        # self.ax.set_ylim(0, None)

    def _setup_cax(self):
        """ Colorbar location setup

        Should be run after creating the figure and before adding data
        """
        # Checks
        self._check_if_closed()
        self._check_if_filled()

        # Make space for colorbar
        div = make_axes_locatable(self.ax)
        self.cax = div.append_axes(
            position="right",
            size="5%",
            pad="5%"
        )

    def add_plot(
        self,
        dataset: xr.Dataset,
        x: Numeric,
        label: str = "",
    ):
        """ Adds data to a plot

        Input:
            dataset:    dataset containing model output
            x:          x-coordinate

        Options:
            label:      label for plot
        """
        # Checks
        self._check_if_closed()
        self._check_if_filled()

        self.x = x

        # Compute spectrum
        wavenumber, freqs, power = fa.spectral_analysis_2d(
            dataset, x=self.x, variable=self.variable, demean=self.demean
        )

        print(freqs)

        # Compute plot limits
        rows = np.any(power > 0.01 * np.max(power), axis=1)
        cols = np.any(power > 0.01 * np.max(power), axis=0)

        # Select cols and rows
        wavenumber = wavenumber[cols]
        freqs = freqs[rows]
        power = power[rows, :][:, cols]

        # Plot data
        self.im = self.ax.pcolormesh(
            self.x_scale * wavenumber,
            self.y_scale * freqs,
            power,
            shading="nearest", # "gouraud",
            cmap=cmo.cm.matter,
            vmin=5.0 * np.min(power)
        )

        # Colorbar
        self.cbar = self.fig.colorbar(self.im, cax=self.cax)
        self.cbar.set_label("Spectral power")

        # Log
        self.data_label = label
        dataset_name: str
        try:
            dataset_name = dataset.attrs["name"]
            if self.data_label == "":
                self.data_label = dataset_name
        except KeyError:
            dataset_name = "'unnamed'"
        print(
            f"# Added data from {dataset_name} for {x=}"
        )
        return self

    def add_dispersion(
        self,
        n: Numeric,
        alpha: Numeric = 1./400.
    ):
        """ Adds the dispersion relation to a plot

        Input:
            n:          maximum number of modes

        Options:
            alpha:      slope for computing dispersion relation
        """
        # Get limits from figure
        xlims = self.ax.get_xlim()
        ylims = self.ax.get_ylim()

        # Compute and plot dispersion relation
        wavenumber = np.linspace(np.min(xlims), np.max(xlims), 100) / self.x_scale

        for i in range(n):
            dispersion = fa.dispersion_relation(wavenumber, n=i, alpha=alpha)
            self.ax.plot(
                self.x_scale * wavenumber,
                self.y_scale * dispersion,
                linewidth=1,
                linestyle="--",
                color="black"
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
        """ Saves the figure

        Input:
            saveloc:    location where the figure should be saved

        Options:
            close:      close figure
        """
        # Checks
        self._check_if_closed()
        savename = f"{saveloc}/spectrum_2d_{plot_spectrum_2d.number:02.0f}" \
            .replace("//", "/")

        # Setup
        self._setup_figure()
        self._setup_plot()

        # Save
        self.fig.savefig(
            savename,
            bbox_inches="tight",
            dpi=FIG_DPI,
            pil_kwargs=FIG_PIL_KWARGS,
        )

        # End
        print(f"# Saved 2d-spectrum figure as {savename}")
        if close:
            self.close()
        return

    def close(self):
        """ Close the figure """
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure 2d-spectrum is closed")
        return


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

    def test_spectrum_2d():
        x = 1e4
        plot_spectrum_2d(variable="wl", demean=True) \
            .add_plot(data_a, x=x) \
            .add_dispersion(n=3) \
            .save(figure_dir)

        plot_spectrum_2d(variable="p", demean=False) \
            .add_plot(data_a, x=x) \
            .save(figure_dir)

    # fmt: on
    test_spectrum_1d()
    test_spectrum_2d()
