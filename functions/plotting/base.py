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
    def __init__(self):
        """Base class for plots containing common methods"""
        self.figsize = FIGSIZE_NORMAL
        self.closed = False
        self.filled = False
