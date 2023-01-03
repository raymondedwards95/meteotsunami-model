""" Functions for contour plots of model output

Main classes:
    plot_contour
"""

import os
import sys

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
from functions.plotting.base import plot_base
import functions.utilities as fu
# fmt: on


class plot_contour(plot_base):
    """Methods to create visualisations of time-slices"""

    number = 0

    def __init__(
        self,
        scale: str = "Mm",
        title: str = None,
    ):
        """Create and setup a figure for time-slices

        Options:
            `scale`:        scale of plots ('m', 'km' or 'Mm')
            `title`:        figure title

        Methods:
            `add_plots`:    add data to the figure
            `save`:         writes the figure to disk
        """
        plot_contour.number += 1
        print("\n# Creating new figure")

        super().__init__()
        self.figsize = FIGSIZE_LONG

        self.title = title
        self.figure_type = "Contours"
        self.figure_num = plot_contour.number
        print(f"# Initiated figure '{self.figure_type} {self.figure_num}'")

        self.x_max = None
        self.x_min = None
        self.y_max = None
        self.y_min = None

        self.unit = None
        self.scale_factor = None
        self.set_scale(scale)

        self._check_if_closed()

    def _setup_figure(self):
        """Figure setup"""
        # Checks
        self._check_if_closed()

        # General
        super()._base_setup_figure()

    def _setup_plot(self):
        """Plot setup"""
        # Checks
        self._check_if_closed()

        # All
        for ax in self.axes.ravel():
            ax.set_xlim(
                fu.none_multiply(self.x_min, 1.0 / self.scale_factor),
                fu.none_multiply(self.x_max, 1.0 / self.scale_factor),
            )
            ax.set_ylim(
                fu.none_multiply(self.y_min, 1.0 / self.scale_factor),
                fu.none_multiply(self.y_max, 1.0 / self.scale_factor),
            )
            # ax.ticklabel_format(scilimits=(-2, 2), useMathText=True)
            ax.grid()

        # Left
        self.fig.supylabel(f"$y$ [{self.unit}]")

        # Bottom
        for j in range(self.axes.shape[1]):
            self.axes[-1, j].set_xlabel(f"$x$ [{self.unit}]")

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
                    f"'{self.figure_type} {self.figure_num}' - "
                    + f"{variable=} should be 'wl', 'u', 'v' or 'p'"
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
        scale: str = None,
    ):
        """Adds a subplot with data

        Input:
            `dataset`:          dataset containing model output
            `variable_list`:    list of names of variables, e.g. 'wl', 'u', 'v' or 'p'
            `t_list`:           list of t-coordinates

        Options:
            `x_min`:            lower limit for x (in meters)
            `x_max`:            upper limit for x (in meters)
            `y_min`:            lower limit for y (in meters)
            `y_max`:            upper limit for y (in meters)
            `scale`:            scale of plots ('m', 'km' or 'Mm')
        """
        # Checks
        self._check_if_closed()
        self._check_if_filled()
        print(f"# Adding plots to '{self.figure_type} {self.figure_num}'...")

        t_list = np.sort(np.array(t_list))
        variable_list = np.char.lower(np.char.strip(np.array(variable_list)))

        if not np.all(np.isin(variable_list, np.array(["wl", "u", "v", "p"]))):
            raise ValueError(
                f"'{self.figure_type} {self.figure_num}' - "
                + f"Variables in `variable_list` should be 'wl', 'u', 'v' or 'p'"
            )

        if scale is not None:
            self.set_scale(scale)

        self.x_max = x_max
        self.x_min = x_min
        self.y_max = y_max
        self.y_min = y_min

        # Create subplots
        var_num = len(variable_list)
        t_num = len(t_list)
        print(f"# Variables:", *variable_list)
        print(f"# t:", *t_list)

        self.fig, self.axes = plt.subplots(
            t_num,
            var_num,
            sharex=True,
            sharey=True,
            squeeze=False,
        )

        # Find plot limits
        var_max = np.zeros(var_num)
        for var_idx, var in enumerate(variable_list):
            var_max[var_idx] = fu.relative_ceil(np.nanmax(np.fabs(dataset[var])))

        # Set colormaps
        cmap_list = []
        norm_list = []
        im_list = []
        for var_idx, var in enumerate(variable_list):
            _var_max = var_max[var_idx]
            _var_min = -1.0 * _var_max
            cmap_list.append(self._select_cmap(var))
            norm_list.append(Normalize(_var_min, _var_max))
            im_list.append(cm.ScalarMappable(norm=norm_list[-1], cmap=cmap_list[-1]))

        # Fill subplots
        for t_idx, t in enumerate(t_list):
            for var_idx, var in enumerate(variable_list):
                print(f"# Adding data for {var:2} and t={t}")
                cont = self.axes[t_idx, var_idx].pcolormesh(
                    dataset["x"] / self.scale_factor,
                    dataset["y"] / self.scale_factor,
                    dataset[var].interp(t=fu.to_timestr(t)),
                    # levels=101,
                    cmap=cmap_list[var_idx],
                    vmin=(-1.0 * var_max[var_idx]),
                    vmax=var_max[var_idx],
                    rasterized=True,
                )

                self.axes[t_idx, var_idx].annotate(
                    f"$t = {t / 3600:0.1f}$h",
                    xy=(0.99, 0.98),
                    xycoords="axes fraction",
                    ha="right",
                    va="top",
                )

        # Add colorbars
        for var_idx, var in enumerate(variable_list):
            cb = self.fig.colorbar(
                im_list[var_idx],
                ax=self.axes[:, var_idx],
                location="top",
                orientation="horizontal",
                label=self._pick_label(var),
                shrink=0.8,
                pad=0.01,
            )
            cb.ax.xaxis.set_label_position("top")
            cb.ax.xaxis.set_ticks_position("bottom")

        # End
        self.filled = True
        return self

    def save(
        self,
        saveloc: str,
        close: bool = True,
    ):
        """Saves the figure

        Input:
            saveloc:    location where the figure should be saved

        Options:
            close:      close figure
        """
        # Checks
        self._check_if_closed()
        savename = f"{saveloc}/contour_{plot_contour.number:02.0f}".replace("//", "/")

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
    data_a = xr.open_dataset(data_a, chunks={"t": "auto", "x": -1, "y": -1})

    # Make figures
    # fmt: off
    def test_contours():
        plot_contour(scale="Mm") \
            .add_plots(
                dataset=data_a,
                variable_list=["wl"],
                t_list=[10 * 3600 * i for i in range(1, 5)],
                y_min=0,
            ) \
            .save(figure_dir)

        plot_contour() \
            .add_plots(
                dataset=data_a,
                variable_list=["u", "v"],
                t_list=[12 * 3600 * i for i in range(1, 4)],
                x_max=4e5,
            ) \
            .save(figure_dir)

        plot_contour() \
            .add_plots(
                dataset=data_a,
                variable_list=["wl", "p"],
                t_list=[10 * 3600],
                x_max=5e5,
                y_min=-5e5,
                y_max=35e5,
                scale="km",
            ) \
            .save(figure_dir)
    # fmt: on

    test_contours()
