""" Scripts to make figures for help with explaining theory """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


### Settings
show_figures = False
current_dir = os.path.dirname(os.path.realpath(__file__))


### Parameters
g = 9.81
a = np.arange(0, 301, 50) * 1e3
a = np.geomspace(1, 301, 50) * 1e3
alpha = 1. / np.logspace(1.5, 3.1)
U = np.arange(0, 60+1)


### Critical Velocity
savename = f"{current_dir}/contour_crit_vel"

plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)
plt.title("Critical Velocity")
Ucr = np.zeros((a.size, alpha.size))
for i in range(a.size):
    for j in range(alpha.size):
        Ucr[i, j] = np.sqrt(g * a[i] * alpha[j] / np.pi)
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
plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
print(f"Saved figure {savename}")


### Wavelength fundamental mode
savename = f"{current_dir}/contour_fundamental_wl"

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
plt.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
print(f"Saved figure {savename}")


### End
if show_figures:
    plt.show()
