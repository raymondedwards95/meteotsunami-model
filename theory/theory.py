""" Scripts to make figures for help with explaining theory """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np


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
plt.savefig(f"{current_dir}/theory_speed_size", bbox_inches="tight")


plt.figure()
for i in range(alpha.size):
    plt.plot(a / 1000., wavelength(u_crit(a, alpha[i]), alpha[i]) / 1000., label=f"$\\alpha = {alpha[i]}$")
plt.plot(a / 1000., 2. * a / 1000., color="black", linestyle="--", linewidth=1, label="$2a$")
plt.legend()
plt.title("Wavelength of Edge Wave Packet")
plt.xlabel("$a$ [km]")
plt.ylabel("$\\lambda$ [km]")
plt.grid()
plt.savefig(f"{current_dir}/theory_wavelength_size", bbox_inches="tight")


plt.show()
