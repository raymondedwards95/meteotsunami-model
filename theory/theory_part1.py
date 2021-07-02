""" Scripts to make figures for help with explaining theory """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    create_map = True
except:
    create_map = False

sns.set_palette(sns.color_palette("muted"))
current_dir = os.path.dirname(os.path.realpath(__file__))


g = 9.81
a = np.logspace(2, 6)
alpha = np.array([1/400, 1/40, 1/4])


def u_crit(a, alpha):
    return np.sqrt(g * a * alpha / np.pi)

def wavelength(U, alpha):
    return 2. * np.pi * U * U / g / alpha


plt.figure()
for i in range(alpha.size):
    plt.semilogx(a / 1000., u_crit(a, alpha[i]), label=f"$\\alpha = {alpha[i]}$")
plt.legend()
plt.title("Critical Storm Speed")
plt.xlabel("$a$ [km]")
plt.ylabel("$U_{crit}$ [m/s]")
plt.grid()
plt.ylim(0, 80)
plt.savefig(f"{current_dir}/line_speed_size", bbox_inches="tight")


plt.figure()
for i in range(alpha.size):
    plt.plot(a / 1000., wavelength(u_crit(a, alpha[i]), alpha[i]) / 1000., label=f"$\\alpha = {alpha[i]}$")
plt.plot(a / 1000., 2. * a / 1000., color="black", linestyle="--", linewidth=1, label="$2a$")
plt.legend()
plt.title("Wavelength of Edge Wave Packet")
plt.xlabel("$a$ [km]")
plt.ylabel("$\\lambda$ [km]")
plt.grid()
plt.savefig(f"{current_dir}/line_wavelength_size", bbox_inches="tight")


if create_map:
    fig = plt.figure()
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
    plt.savefig(f"{current_dir}/map", bbox_inches="tight")


plt.show()
