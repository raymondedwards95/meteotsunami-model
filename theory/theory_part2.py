""" Scripts to make figures for help with explaining the theory parts """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.theory as ft
# fmt: on


# Critical Velocity
def theory_figure_speed_vs_size_alpha(
    a: npt.ArrayLike, alpha: npt.ArrayLike, savename: str = None
) -> None:
    """Plots the relation between critical wave velocity and pressure disturbance size and the slope of a flat sloped bottom topography

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
    Ucr = ft.critical_velocity_sloped(grid_a, grid_alpha)

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Critical Velocity", va="top", ha="left", x=0.01)

    cl = plt.contour(
        alpha,
        a / 1e3,
        Ucr,
        levels=np.arange(U.min(), U.max() + 10, 10),
        colors="black",
        alpha=0.8,
    )
    plt.clabel(cl, fmt="%2.0f")
    plt.contourf(alpha, a / 1e3, Ucr, levels=U, cmap=cmo.cm.speed)
    cb = plt.colorbar()
    cb.ax.set_title("$U_{cr}$ [m/s]")
    plt.xlabel("$\\alpha$ [-]")
    plt.ylabel("$a$ [km]")
    plt.grid()
    plt.xscale("log")
    for _angle in [1 / 50, 1 / 400, 1 / 200, 1 / 1000]:
        plt.text(
            _angle,
            0.01,
            f"1/{1/_angle:0.0f}",
            color="red",
            transform=plt.gca().get_xaxis_transform(),
        )
        plt.axvline(_angle, color="red")
    fig.get_layout_engine().execute(fig)
    fig.savefig(
        savename,
        bbox_inches="tight",
        dpi=FIG_DPI,
        pil_kwargs=FIG_PIL_KWARGS,
    )

    print(f"Saved figure {savename}")
    return


# Wavelength fundamental mode
def theory_figure_wavelength_vs_alpha_speed(
    alpha: npt.ArrayLike, velocity: npt.ArrayLike, savename: str = None
) -> None:
    """Plots the relation between wavelength and bottom slope and pressure disturbance velocity

    Input:
        `alpha`:    array with bottom slope
        `U`:        array with velocity of pressure disturbance

    Options:
        `savename`  figure name
    """
    # Arguments
    if savename is None:
        savename = f"{figure_dir}/contour_fundamental_wl"

    # Computations
    grid_velocity, grid_alpha = np.meshgrid(velocity, alpha, sparse=True)
    wavelength = ft.fundamental_wavelength_sloped(grid_velocity, grid_alpha).T

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Wavelength Fundamental Mode", va="top", ha="left", x=0.01)

    cl = plt.contour(
        alpha,
        velocity,
        wavelength / 1e3,
        levels=np.arange(0, 4e4 + 10, 1e4) / 1e3,
        colors="black",
        alpha=0.8,
    )
    plt.clabel(cl, fmt="%2.0f")
    plt.contourf(
        alpha,
        velocity,
        wavelength / 1e3,
        levels=np.linspace(0, 40000, 51) / 1e3,
        cmap=cmo.cm.amp,
    )
    cb = plt.colorbar()
    cb.ax.set_title("$\\lambda_0$ [km]")
    plt.xlabel("$\\alpha$ [-]")
    plt.ylabel("$U$ [m/s]")
    plt.grid()
    plt.xscale("log")
    for _angle in [1 / 50, 1 / 400, 1 / 200, 1 / 1000]:
        plt.text(
            _angle,
            0.01,
            f"1/{1/_angle:0.0f}",
            color="blue",
            transform=plt.gca().get_xaxis_transform(),
        )
        plt.axvline(_angle, color="blue")
    fig.get_layout_engine().execute(fig)
    plt.savefig(
        savename,
        bbox_inches="tight",
        dpi=FIG_DPI,
        pil_kwargs=FIG_PIL_KWARGS,
    )

    print(f"Saved figure {savename}")
    return


if __name__ == "__main__":
    # Paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    figure_dir = f"{script_dir}/figures"
    os.makedirs(figure_dir, exist_ok=True)

    # Parameters
    # a = np.arange(0, 301, 50) * 1e3
    a = np.geomspace(1, 301, 100) * 1e3
    # alpha = 1. / np.logspace(1.5, 3.1, 100)
    alpha = np.geomspace(1e-4, 1e-1, 100)
    U = np.arange(0, 60 + 0.1, 0.5)

    # Figures
    theory_figure_speed_vs_size_alpha(a, alpha)
    theory_figure_wavelength_vs_alpha_speed(alpha, U)
