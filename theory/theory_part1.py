""" Scripts to make figures for help with explaining the theory parts """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scipy.constants

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    create_map = True
except:
    print("Cartopy is not installed!")
    create_map = False

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


### Functions
def u_crit(a, alpha):
    """ Computes critical velocity of pressure perturbation for given size `a` and bottom slope `alpha` """
    return np.sqrt(g * a * alpha / np.pi)

def wavelength(U, alpha):
    """ Computes wavelength for given velocity `U` and bottom slope `alpha` """
    return 2. * np.pi * U * U / g / alpha


### Speed
def theory_figure_speed_vs_size(a: npt.ArrayLike, savename: str=None) -> None:
    if savename is None:
        savename = f"{figure_dir}/line_speed_size"

    plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)
    for i in range(alpha.size):
        plt.semilogx(a / 1000., u_crit(a, alpha[i]), label=f"$\\alpha = {alpha[i]}$")
    plt.legend()
    plt.title("Critical Storm Speed")
    plt.xlabel("$a$ [km]")
    plt.ylabel("$U_{crit}$ [m/s]")
    plt.grid()
    plt.ylim(0, 80)
    plt.xlim(1e-1, 1e3)
    plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)

    print(f"Saved figure {savename}")
    return


### Wavelength
def theory_figure_wavelength_vs_size(a: npt.ArrayLike, savename:str=None) -> None:
    if savename is None:
        savename = f"{figure_dir}/line_wavelength_size"

    plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)
    for i in range(alpha.size):
        plt.plot(a / 1000., wavelength(u_crit(a, alpha[i]), alpha[i]) / 1000., label=f"$\\alpha = {alpha[i]}$")
    plt.plot(a / 1000., 2. * a / 1000., color="black", linestyle="--", linewidth=1, label="$2a$")
    plt.legend()
    plt.title("Wavelength of Edge Wave Packet")
    plt.xlabel("$a$ [km]")
    plt.ylabel("$\\lambda$ [km]")
    plt.xlim(0, 1000)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, 1001, 100))
    plt.yticks(np.arange(0, 2001, 200))
    plt.grid()
    plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)

    print(f"Saved figure {savename}")
    return


### Map
def theory_figure_map(savename: str=None) -> None:
    if savename is None:
        savename = f"{figure_dir}/map"

    fig = plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)
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

    a = np.logspace(2, 6)
    alpha = np.array([1/400, 1/40, 1/4])


    # Figures
    theory_figure_speed_vs_size(a)
    theory_figure_wavelength_vs_size(a)
    if create_map:
        theory_figure_map()


    # End
    if show_figures:
        plt.show()
