""" Functions for contour plots of model output

Main functions:
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
# fmt: on


class plot_contour():
    """ Methods to create visualisations of time-slices """
    number = 0

    def __init__(self):
        """ Create and setup a figure for time-slices
        """
        raise NotImplementedError

    def _check_if_closed(self):
        """ Raises an error if the figure is supposed to be closed """
        if self.closed:
            raise TypeError(
                f"Figure is cleared and closed: it can not be edited")

    def _setup_figure(self):
        """ Figure setup """
        # Checks
        self._check_if_closed()
        raise NotImplementedError

    def _setup_plot(self):
        """ Plot setup """
        # Checks
        self._check_if_closed()
        raise NotImplementedError

    def add_plot(
        self,
        dataset: xr.Dataset,
        variable: str,
        t: Numeric,
        ):
        """ Adds a subplot with data

        Input:
            dataset:    dataset containing model output
            variable:   name of variable, e.g. 'wl', 'u', 'v' or 'p'
            x:          x-coordinate
        """
        raise NotImplementedError

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
