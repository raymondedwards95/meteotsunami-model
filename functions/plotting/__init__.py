""" Functions for visualising model output """
import os
import sys

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
from functions.plotting.spectrum import plot_spectrum_1d
# fmt: on
