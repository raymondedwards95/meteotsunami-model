""" Functions for contour plots of model output

Main functions:
"""

import os
import sys

import cmocean as cmo
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import xarray as xr
from matplotlib.colors import Normalize

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
import functions.utilities as fu
# fmt: on


class plot_contour():
    """ Methods to create visualisations of time-slices """
    number = 0

    def __init__(self):
        """ Create and setup a figure for time-slices
        """
        plot_contour.number += 1

        self.figsize = FIGSIZE_LONG
        self.closed = False
        self.filled = False

        self.x_max = None
        self.x_min = None
        self.y_max = None
        self.y_min = None

        print(
            f"\n# Initiated figure for contour-levels"
        )

    def _check_if_closed(self):
        """ Raises an error if the figure is supposed to be closed """
        if self.closed:
            raise TypeError(
                f"Figure is cleared and closed: it can not be edited"
            )

    def _check_if_filled(self):
        """ Raises an error if the figure is filled with data already """
        if self.filled:
            raise TypeError(
                f"Figure contains data already: it cannot be edited"
            )

    def _setup_figure(self):
        """ Figure setup """
        # Checks
        self._check_if_closed()

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        # self.fig.set_tight_layout(True)

        # Figure specific
        self.fig.suptitle(
            f"Contours", va="top", ha="left", x=0.01
        )

    def _setup_plot(self):
        """ Plot setup """
        # Checks
        self._check_if_closed()
        print(self.axes.shape)

        # All
        for ax in self.axes.ravel():
            ax.set_xlim(self.x_min, self.x_max)
            ax.set_ylim(self.y_min, self.y_max)

        # Left
        for i in range(self.axes.shape[0]):
            self.axes[i, 0].set_ylabel("$y$ [km]")

        # Bottom
        for j in range(self.axes.shape[1]):
            self.axes[-1, j].set_xlabel("$x$ [km]")

    def _pick_cmap(self, variable: str):
        match variable.lower().strip():
            case "wl":
                return cmo.cm.balance
            case ("u" | "v"):
                return cmo.cm.delta
            case "p":
                return cmo.cm.curl
            case _:
                raise ValueError(
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

    def _pick_label(self, variable: str):
        match variable.lower().strip():
            case "wl":
                return "Water level [m]"
            case "u":
                return "Cross shore water velocity [m/s]"
            case "v":
                return "Along shore water velocity [m/s]"
            case "p":
                return "Surface air pressure [Pa]"
            case _:
                raise ValueError(
                    f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

    def add_plots(
        self,
        dataset: xr.Dataset,
        variable_list: npt.ArrayLike,
        t_list: npt.ArrayLike,
        x_min: Numeric = None,
        x_max: Numeric = None,
        y_min: Numeric = None,
        y_max: Numeric = None,
    ):
        """ Adds a subplot with data

        Input:
            `dataset`:          dataset containing model output
            `variable_list`:    list of names of variables, e.g. 'wl', 'u', 'v' or 'p'
            `t_list`:           list of t-coordinates

        Options:
            `x_max`:            upper limit for x (in kilometers)
            `x_min`:            lower limit for x (in kilometers)
            `y_max`:            upper limit for y (in kilometers)
            `y_min`:            lower limit for y (in kilometers)
        """
        # Checks
        self._check_if_closed()
        self._check_if_filled()
        print(f"# Adding plots to figure...")

        t_list = np.sort(np.array(t_list))
        variable_list = np.char.lower(np.char.strip(np.array(variable_list)))

        if not np.all(np.isin(variable_list, np.array(["wl", "u", "v", "p"]))):
            raise ValueError(
                f"Variables in `variable_list` should be 'wl', 'u', 'v' or 'p'"
            )

        self.x_max = x_max
        self.x_min = x_min
        self.y_max = y_max
        self.y_min = y_min

        # Create subplots
        var_num = len(variable_list)
        t_num = len(t_list)
        print(f"# Variables:", *variable_list)
        print(f"# t:        ", *t_list)

        self.fig, self.axes = plt.subplots(
            t_num,
            var_num,
            sharex=True,
            sharey=True,
            squeeze=False,
            constrained_layout=True,
        )

        # Find plot limits
        var_max = np.zeros(var_num)
        for var_idx, var in enumerate(variable_list):
            var_max[var_idx] = np.nanmax(np.fabs(dataset[var]))

        # Set colormaps
        cmap_list = []
        norm_list = []
        im_list = []
        for var_idx, var in enumerate(variable_list):
            _var_max = var_max[var_idx]
            _var_min = -1. * _var_max
            cmap_list.append(self._pick_cmap(var))
            norm_list.append(Normalize(_var_min, _var_max))
            im_list.append(cm.ScalarMappable(
                norm=norm_list[-1], cmap=cmap_list[-1]))

        # Fill subplots
        for t_idx, t in enumerate(t_list):
            for var_idx, var in enumerate(variable_list):
                print(f"# Adding data for {var:2} and t={t}")
                self.axes[t_idx, var_idx].contourf(
                    dataset["x"] / 1000.,
                    dataset["y"] / 1000.,
                    dataset[var].interp(t=fu.to_timestr(t)),
                    levels=101,
                    cmap=cmap_list[var_idx],
                    vmin=(-1. * var_max[var_idx]),
                    vmax=var_max[var_idx],
                )
                self.axes[t_idx, var_idx].annotate(
                    f"$t = {t / 3600:0.1f}$h",
                    xy=(0.99, 0.98),
                    xycoords="axes fraction",
                    ha="right",
                    va="top"
                )

        # Add colorbars
        for var_idx, var in enumerate(variable_list):
            self.fig.colorbar(
                im_list[var_idx],
                ax=self.axes[:, var_idx],
                location="top",
                label=self._pick_label(var),
                aspect=50
            )

        # End
        self.filled = True
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
        savename = f"{saveloc}/contour_{plot_contour.number:02.0f}" \
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
        print(f"# Saved contour figure as {savename}")
        if close:
            self.close()
        return

    def close(self):
        """ Close the figure """
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure contour is closed")
        return


if __name__ == "__main__":
    # Define paths
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    figure_dir = f"{script_dir}/tests/figures/"
    os.makedirs(figure_dir, exist_ok=True)

    # Get data
    main_dir = os.path.dirname(script_dir)
    data_a = f"{main_dir}/reproduction-an-2012/output/data_repr_00.nc"

    data_a = xr.open_dataset(data_a, chunks={"t": "auto", "x": -1, "y": -1})

    # Make figures
    def test_contours():
        plot_contour() \
            .add_plots(
                dataset=data_a,
                variable_list=["wl"],
                t_list=[10 * 3600 * i for i in range(1, 5)],
        ) \
            .save(figure_dir)

        plot_contour() \
            .add_plots(
                dataset=data_a,
                variable_list=["u", "v"],
                t_list=[10 * 3600 * i for i in range(1, 5)],
                x_max=400,
        ) \
            .save(figure_dir)

        plot_contour() \
            .add_plots(
                dataset=data_a,
                variable_list=["wl", "p"],
                t_list=[10 * 3600],
                x_max=500,
                y_min=-500,
                y_max=3500,
        ) \
            .save(figure_dir)

    test_contours()
