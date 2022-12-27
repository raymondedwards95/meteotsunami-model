""" Common functions for plots of model output

Main classes:
    plot_base
"""

import os
import sys

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

    def _check_if_closed(self):
        """Raises an error if the figure is supposed to be closed"""
        if self.closed:
            raise TypeError(
                f"Figure '{self.figure_type}' is cleared and closed: it can not be edited"
            )

    def _check_if_filled(self):
        """Raises an error if the figure is filled with data already"""
        if self.filled:
            raise TypeError(
                f"Figure '{self.figure_type}' contains data already: it cannot be edited"
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
                raise ValueError(f"{variable=} should be 'wl', 'u', 'v' or 'p'")

        return self
