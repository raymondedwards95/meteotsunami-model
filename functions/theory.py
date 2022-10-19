""" Theoretical relationships

Main functions:
    critical_velocity_sloped
"""

import os
import sys

import numpy as np
import scipy.constants

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.utilities as fu
# fmt: on


g = scipy.constants.g
assert np.abs(g - 9.81) < 1e-2


@np.vectorize
def critical_velocity_sloped(a: Numeric, alpha: Numeric) -> Numeric:
    """ Computes the critical wave velocity

    Input:
        `a`:        pressure disturbance size
        `alpha`:    bottom slope

    Output:
        `Ucr`:      critical wave velocity
    """
    return np.sqrt(g * a * alpha / np.pi)


@np.vectorize
def fundamental_wavelength_sloped(velocity: Numeric, alpha: Numeric) -> Numeric:
    """ Computes the wavelength of the fundamental mode

    Input:
        `velocity`:     wave velocity
        `alpha`:        bottom slope

    Output:
        `wavelength`:   critical wave velocity
    """
    return 2. * np.pi * velocity / g / alpha
