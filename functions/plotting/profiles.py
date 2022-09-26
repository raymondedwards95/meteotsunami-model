""" Functions for visualising wave profiles from model output

Main classes:
    plot_alongshore
    plot_crossshore
"""

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
import functions.analysis as fa
# fmt: on


class plot_alongshore():
    raise NotImplementedError


class plot_crossshore():
    raise NotImplementedError
