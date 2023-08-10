""" Functions for visualising model output """

import os
import sys

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
from functions.plotting.alongshore import plot_alongshore
from functions.plotting.base import plot_base
from functions.plotting.contour import plot_contour
from functions.plotting.crossshore import plot_crossshore
from functions.plotting.growth import plot_growth
from functions.plotting.spectrum1d import plot_spectrum_1d
from functions.plotting.spectrum2d import plot_spectrum_2d
from functions.plotting.timeseries import plot_parametric, plot_timeseries
# fmt: on


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", f"{PATH_PLOTTING}/base.py"])

    subprocess.run(["python", f"{PATH_PLOTTING}/alongshore.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/contour.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/crossshore.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/growth.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/spectrum1d.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/spectrum2d.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/timeseries.py"])
