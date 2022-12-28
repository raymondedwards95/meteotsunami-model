""" Functions for visualising model output """

import os
import sys

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
from functions.plotting.base import plot_base
from functions.plotting.contour import plot_contour
from functions.plotting.profiles import plot_alongshore, plot_crossshore
from functions.plotting.spectrum1d import plot_spectrum_1d
from functions.plotting.spectrum2d import plot_spectrum_2d
from functions.plotting.timeseries import plot_timeseries
# fmt: on


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", f"{PATH_PLOTTING}/base.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/contour.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/profiles.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/spectrum1d.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/spectrum2d.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/timeseries.py"])
