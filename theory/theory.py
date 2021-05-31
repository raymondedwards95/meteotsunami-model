""" Scripts to make figures for help with explaining theory """


import matplotlib.pyplot as plt
import numpy as np

g = 9.81
a = np.arange(0, 15) * 1e5
alpha = np.array([1/400, 1/40])


def u_crit(a, alpha):
    return np.sqrt(g * a * alpha / np.pi)

def wavelength(U, alpha):
    return 2. * np.pi * U * U / g / alpha


plt.figure()
for i in range(alpha.size):
    plt.plot(a / 1000., u_crit(a, alpha[i]), label=f"$\\alpha = {alpha[i]}$")
plt.legend()
plt.title("Critical Storm Speed")
plt.xlabel("$a$ [km]")
plt.ylabel("$U_{crit}$ [m/s]")


plt.figure()
for i in range(alpha.size):
    plt.plot(a / 1000., wavelength(u_crit(a, alpha[i]), alpha[i]) / 1000., label=f"$\\alpha = {alpha[i]}$")
plt.plot(a / 1000., 2. * a / 1000., color="black", linestyle="--", linewidth=1, label="$2a$")
plt.legend()
plt.title("Wavelength of Edge Wave Packet")
plt.xlabel("$a$ [km]")
plt.ylabel("$\\lambda$ [km]")


plt.show()
