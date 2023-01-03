""" Common functions for plots of model output

Main classes:
    plot_base
"""

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
# fmt: on


class plot_base:
    """Base class for plots containing common methods"""

    def __init__(self):
        self.figsize = FIGSIZE_NORMAL
        self.closed = False
        self.filled = False
        self.title = None
        self.figure_type = "Base"
        self.figure_num = -1
        self.fig = plt.figure()

    def _check_if_closed(self):
        """Raises an error if the figure is supposed to be closed"""
        if self.closed:
            raise TypeError(
                f"Figure '{self.figure_type} {self.figure_num}' is cleared and closed: "
                + f"it can not be edited"
            )

    def _check_if_filled(self):
        """Raises an error if the figure is filled with data already"""
        if self.filled:
            raise TypeError(
                f"Figure '{self.figure_type} {self.figure_num}' contains data already: "
                + f"it cannot be edited"
            )

    def _base_setup_figure(self):
        """Base settings for figure"""
        # Checks
        self._check_if_closed()

        # General
        self.fig.set_size_inches(self.figsize)
        self.fig.set_dpi(FIG_DPI)
        self.fig.set_layout_engine("compressed")

        # Title
        if self.title is None:
            self.title = self.figure_type
        self.fig.suptitle(self.title, va="top", ha="left", x=0.01)

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
                raise ValueError(
                    f"'{self.figure_type} {self.figure_num}' - "
                    + f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

        print(
            f"# Variable in '{self.figure_type} {self.figure_num}'",
            f"is '{self.variable}' with units '{self.variable_unit}'",
        )

        return self

    def _select_cmap(self, variable: str):
        match variable.lower().strip():
            case "wl":
                return cmo.cm.balance
            case ("u" | "v"):
                return cmo.cm.delta
            case "p":
                return cmo.cm.curl
            case _:
                raise ValueError(
                    f"'{self.figure_type} {self.figure_num}' - "
                    + f"{variable=} should be 'wl', 'u', 'v' or 'p'"
                )

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
                    f"Scale should be either 'm', 'km', or 'Mm' "
                    + f"instead of '{scale}'"
                )

        print(
            f"# Set scale of '{self.figure_type} {self.figure_num}'",
            f"to {self.unit} (=1/{self.scale_factor})",
        )

        return self

    def close(self):
        """Close the figure"""
        self.closed = True
        self.fig.clear()
        plt.close(self.fig)
        print(f"# Figure '{self.figure_type} {self.figure_num}' is closed")
        return


if __name__ == "__main__":
    figure = plot_base()
    figure._check_if_closed()
    for unit in ["m", "km", "Mm"]:
        figure.set_scale(unit)
    figure._check_if_filled()
    for var in ["p", "u", "v", "wl"]:
        figure._set_variable(var)
    figure._base_setup_figure()
    figure.close()
