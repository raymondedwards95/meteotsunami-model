""" Scripts to make figures for help with explaining the theory parts """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scipy.constants

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
from functions.theory import critical_velocity_sloped
# fmt: on


### Critical Velocity
def theory_figure_speed_vs_size_alpha(a: npt.ArrayLike, alpha: npt.ArrayLike, savename: str=None) -> None:
    """ Plots the relation between critical wave velocity and pressure disturbance size and the slope of a flat sloped bottom topography

    Input:
        `a`:        array with pressure disturbance size
        `alpha`:    array with bottom slope

    Options:
        `savename`  figure name
    """
    # Arguments
    if savename is None:
        savename = f"{figure_dir}/contour_crit_vel"

    # Computations
    grid_a, grid_alpha = np.meshgrid(a, alpha, sparse=True)
    Ucr = critical_velocity_sloped(grid_a, grid_alpha)

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Critical Velocity", va="top", ha="left", x=0.01)

    cl = plt.contour(alpha, a/1e3, Ucr, levels=np.arange(U.min(), U.max()+10, 10), colors="black")
    plt.clabel(cl, fmt="%2.0f")
    plt.contourf(alpha, a/1e3, Ucr, levels=U, cmap=cmo.cm.speed)
    cb = plt.colorbar()
    cb.ax.set_title("$U_{cr}$ [m/s]")
    plt.xlabel("$\\alpha$ [-]")
    plt.ylabel("$a$ [km]")
    plt.grid()
    plt.xscale("log")
    for _angle in [1/50, 1/200, 1/1000]:
        plt.text(_angle, 0.01, f"{_angle}", color="red", transform=plt.gca().get_xaxis_transform())
        plt.axvline(_angle, color="red")
    fig.get_layout_engine().execute(fig)
    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)

    print(f"Saved figure {savename}")
    return


### Wavelength fundamental mode
def theory_figure_wavelength_vs_alpha_speed(alpha: npt.ArrayLike, U: npt.ArrayLike, savename: str=None) -> None:
    """ Plots the relation between wavelength and bottom slope and pressure disturbance velocity

    Input:
        `alpha`:    array with bottom slope
        `U`:        array with velocity of pressure disturbance

    Options:
        `savename`  figure name
    """
    if savename is None:
        savename = f"{figure_dir}/contour_fundamental_wl"

    plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)
    plt.title("Wavelength fundamental mode")
    labda = np.zeros((U.size, alpha.size))
    for i in range(U.size):
        for j in range(alpha.size):
            labda[i, j] = 2. * np.pi * U[i] / g / alpha[j]
    cl = plt.contour(alpha, U, labda/1e3, levels=np.arange(0, 4e4+10, 1e4)/1e3, colors="black")
    plt.clabel(cl, fmt="%2.0f")
    plt.contourf(alpha, U, labda/1e3, levels=np.linspace(0, 40000, 51)/1e3, cmap=cmo.cm.amp)
    cb = plt.colorbar()
    cb.ax.set_title("$\\lambda_0$ [km]")
    plt.xlabel("$\\alpha$ [-]")
    plt.ylabel("$U$ [m/s]")
    plt.grid()
    plt.xscale("log")
    for _angle in [1/50, 1/200, 1/1000]:
        plt.text(_angle, 0.01, f"{_angle}", color="orange", transform=plt.gca().get_xaxis_transform())
        plt.axvline(_angle, color="orange")
    plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)

    print(f"Saved figure {savename}")
    return


if __name__ == "__main__":
    # Settings
    show_figures = False
    current_dir = os.path.dirname(os.path.realpath(__file__))
    figure_dir = f"{current_dir}/figures"
    os.makedirs(figure_dir, exist_ok=True)

    # Constants
    g = scipy.constants.g
    assert np.abs(g - 9.81) < 1e-2

    # Parameters
    # a = np.arange(0, 301, 50) * 1e3
    a = np.geomspace(1, 301, 50) * 1e3
    alpha = 1. / np.logspace(1.5, 3.1)
    U = np.arange(0, 60+1)

    # Figures
    theory_figure_speed_vs_size_alpha(a, alpha)
    theory_figure_wavelength_vs_alpha_speed(alpha, U)

    # End
    if show_figures:
        plt.show()
