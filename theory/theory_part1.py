""" Scripts to make figures for help with explaining the theory parts """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    create_map = True
except:
    print("Cartopy is not installed!")
    create_map = False

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.theory as ft
# fmt: on


### Speed
def theory_figure_speed_vs_size(a: npt.ArrayLike, alpha: npt.ArrayLike, savename: str = None) -> None:
    """ Plots the relation between critical velocity and pressure disturbance size

    Input:
        `a`:        array with pressure disturbance size
        `alpha`:    array with bottom slope

    Options:
        `savename`  figure name
    """
    # Arguments
    if savename is None:
        savename = f"{figure_dir}/line_speed_size"

    # Computations
    grid_a, grid_alpha = np.meshgrid(a, alpha, sparse=True)
    u_crit = ft.critical_velocity_sloped(grid_a, grid_alpha)

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Critical Storm Speed", va="top", ha="left", x=0.01)

    for i in range(alpha.size):
        plt.semilogx(a / 1000., u_crit[i,:], label=f"$\\alpha = 1/{1/alpha[i]:.0f} \\approx {alpha[i]}$")
        plt.fill_between(a / 1000., u_crit[i,:], alpha=0.1)
    plt.legend()
    plt.xlabel("$a$ [km]")
    plt.ylabel("$U_{crit}$ [m/s]")
    plt.grid()
    plt.ylim(0, 100)
    plt.xlim(1e0, 1e3)
    fig.get_layout_engine().execute(fig)
    plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)

    print(f"Saved figure {savename}")
    return


### Wavelength
def theory_figure_wavelength_vs_size(a: npt.ArrayLike, alpha: npt.ArrayLike, savename: str=None) -> None:
    """ Plots the relation between wavelength and pressure disturbance size

    Input:
        `a`:        array with pressure disturbance size
        `alpha`:    array with bottom slope

    Options:
        `savename`  figure name
    """
    # Arguments
    if savename is None:
        savename = f"{figure_dir}/line_wavelength_size"

    # Computations
    grid_a, grid_alpha = np.meshgrid(a, alpha, sparse=True)
    u_crit = ft.critical_velocity_sloped(grid_a, grid_alpha)
    wavelength = ft.fundamental_wavelength_sloped(u_crit, grid_alpha)

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Wavelength of Edge Wave Packet", va="top", ha="left", x=0.01)

    for i in range(alpha.size):
        plt.plot(a / 1000., wavelength[i,:] / 1000., label=f"$\\alpha = 1/{1/alpha[i]:.0f} \\approx {alpha[i]}$")
        plt.fill_between(a / 1000., wavelength[i,:] / 1000., alpha=0.1)
    plt.plot(a / 1000., 2. * a / 1000., color="black", linestyle="--", linewidth=1, label="$2a$")
    plt.legend()
    plt.xlabel("$a$ [km]")
    plt.ylabel("$\\lambda$ [km]")
    plt.xticks(np.arange(0, 301, 100))
    plt.yticks(np.arange(0, 601, 200))
    plt.xlim(0, 300)
    plt.ylim(0, 700)
    plt.grid()
    fig.get_layout_engine().execute(fig)
    plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)

    print(f"Saved figure {savename}")
    return


### Map
def theory_figure_map(savename: str=None) -> None:
    """ Create figure of the area of interest

    Options:
        savename    figure name
    """
    # Arguments
    if savename is None:
        savename = f"{figure_dir}/map"

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Area of Interest", va="top", ha="left", x=0.01)

    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    ax.coastlines("50m")
    ax.add_feature(
        cfeature.NaturalEarthFeature(
            "physical", "land", "50m",
            edgecolor="black", facecolor="gray"
        )
    )
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    ax.set_xlim([-10, 10])
    ax.set_ylim([45, 60])
    # plt.xlabel("Longitude")
    # plt.ylabel("Latitude")
    fig.get_layout_engine().execute(fig)
    plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)

    print(f"Saved figure {savename}")
    return


if __name__ == "__main__":
    # Settings
    show_figures = False
    current_dir = os.path.dirname(os.path.realpath(__file__))
    figure_dir = f"{current_dir}/figures"
    os.makedirs(figure_dir, exist_ok=True)

    # Parameters
    a = np.logspace(2, 6)
    alpha = np.array([1/400, 1/40, 1/4])

    # Figures
    theory_figure_speed_vs_size(a, alpha)
    theory_figure_wavelength_vs_size(a, alpha)
    if create_map:
        theory_figure_map()

    # End
    if show_figures:
        plt.show()
