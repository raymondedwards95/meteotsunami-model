""" Scripts to make figures for help with explaining the theory parts """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Critical Velocity", va="top", ha="left", x=0.01)

    # Plots
    cl = ax.contour(
        alpha,
        a / 1e3,
        Ucr,
        levels=np.arange(0, 60 + 0.1, 10)[:-1],
        colors="black",
        alpha=0.8,
    )
    ax.clabel(cl, fmt="%2.0f")

    cf = ax.contourf(
        alpha,
        a / 1e3,
        Ucr,
        levels=np.arange(0, 60 + 0.1, 1),
        cmap=cmo.cm.speed,
    )

    # Colorbar
    div = make_axes_locatable(ax)
    cax = div.append_axes(position="right",size="5%",pad="5%")
    cb = fig.colorbar(cf, cax=cax)
    cb.set_label("$U_{cr}$ [m/s]")
    cb.set_ticks(np.arange(0, 60 + 1, 10))
    cb.set_ticks(np.arange(0, 60 + 1, 5), minor=True)

    # Other
    ax.set_xlabel("$\\alpha$ [-]")
    ax.set_ylabel("$a$ [km]")
    ax.grid()
    ax.set_xscale("log")

    # Lines of constant alpha
    for _angle in [1 / 50, 1 / 400, 1 / 200, 1 / 800]:
        ax.text(
            _angle,
            0.01,
            f"1/{1/_angle:0.0f}",
            color="red",
            transform=ax.get_xaxis_transform(),
        )
        ax.axvline(_angle, color="red")

    # Save figure
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
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Wavelength Fundamental Mode", va="top", ha="left", x=0.01)

    # Plots
    cl = ax.contour(
        alpha,
        velocity,
        wavelength / 1e3,
        levels=np.geomspace(5000, 5e5, 3) / 1e3,
        colors="black",
        alpha=0.8,
    )
    ax.clabel(cl, fmt="%2.0f")

    cf = ax.contourf(
        alpha,
        velocity,
        wavelength / 1e3,
        levels=np.arange(0, 7e5 + 1, 1e4) / 1e3,
        cmap=cmo.cm.amp,
        # extend="max",
    )

    # Colorbar
    div = make_axes_locatable(ax)
    cax = div.append_axes(position="right",size="5%",pad="5%")
    cb = fig.colorbar(cf, cax=cax)
    cb.set_label("$\\lambda_0$ [km]")
    cb.set_ticks(np.arange(0, 700 + 1, 100))
    cb.set_ticks(np.arange(0, 500 + 1, 25), minor=True)

    # Other
    ax.set_xlabel("$\\alpha$ [-]")
    ax.set_ylabel("$U$ [m/s]")
    ax.grid()
    ax.set_xscale("log")

    # Lines of constant alpha
    for _angle in [1 / 50, 1 / 400, 1 / 200, 1 / 800]:
        ax.text(
            _angle,
            0.01,
            f"1/{1/_angle:0.0f}",
            color="blue",
            transform=ax.get_xaxis_transform(),
        )
        ax.axvline(_angle, color="blue")

    # Save figure
    fig.get_layout_engine().execute(fig)
    fig.savefig(
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
