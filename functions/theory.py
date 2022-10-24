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
    """Computes the critical wave velocity

    Input:
        `a`:        pressure disturbance size
        `alpha`:    bottom slope

    Output:
        `Ucr`:      critical wave velocity
    """
    return np.sqrt(g * a * alpha / np.pi)


@np.vectorize
def fundamental_wavelength_sloped(velocity: Numeric, alpha: Numeric) -> Numeric:
    """Computes the wavelength of the fundamental mode

    Input:
        `velocity`:     wave velocity
        `alpha`:        bottom slope

    Output:
        `wavelength`:   wavelength of fundamental mode
    """
    return 2.0 * np.pi * np.power(velocity, 2.0) / g / alpha


if __name__ == "__main__":
    a = np.geomspace(1, 1e6, 50)
    alpha = np.geomspace(1e-1, 1e-4, 100)
    vel = np.linspace(0, 100, 30)

    def test_critical_velocity_sloped():
        a_mesh, alpha_mesh = np.meshgrid(a, alpha, sparse=True)

        print(critical_velocity_sloped(a_mesh, alpha_mesh).shape)
        print(critical_velocity_sloped(a.mean(), alpha).shape)
        print(critical_velocity_sloped(a, alpha.mean()).shape)
        print(critical_velocity_sloped(a.mean(), alpha.mean()).shape)

    def test_fundamental_wavelength_sloped():
        vel_mesh, alpha_mesh = np.meshgrid(vel, alpha, sparse=True)

        print(fundamental_wavelength_sloped(vel_mesh, alpha_mesh).shape)
        print(fundamental_wavelength_sloped(vel.mean(), alpha).shape)
        print(fundamental_wavelength_sloped(vel, alpha.mean()).shape)
        print(fundamental_wavelength_sloped(vel.mean(), alpha.mean()).shape)

    test_critical_velocity_sloped()
    test_fundamental_wavelength_sloped()
