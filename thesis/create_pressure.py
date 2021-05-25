""" Creates a file with atmospheric pressure """

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp


### Parameters
# generic
cases = [0]
num_cases = len(cases)

# pressure distribution
t0 = 10000.  # default: 10000 s
U_list = [50.]  # default: 50 m/s
a_list = [200000.]  # default: 200 km
p0_list = [2000.]  # default: 2000 Pa

# cross shore (meters)
x_min = 0.  # default: 0 km
x_max = 1e6  # default: 1000 km
x_step = 10e3  # default: 10 km

# along shore (meters)
y_min = -1e7  # default: -10000 km
y_max = 1e7  # default: 10000 km
y_step = x_step  # default: 10 km

# time (seconds)
t_min = 0  # default: 0
t_max = 70 * 3600.  # default: 70 hours
t_step = 3600. / 4.  # default: 1 hour

# check parameters
assert len(U_list) == num_cases
assert len(a_list) == num_cases
assert len(p0_list) == num_cases

print("\nPressure fields for the following cases are computed (case, U, a, p0)")
for i in range(num_cases):
    print(cases[i], U_list[i], a_list[i], p0_list[i])


### Grid
x_num = int((x_max - x_min) / x_step + 1)
y_num = int((y_max - y_min) / y_step + 1)

x = np.linspace(x_min, x_max, x_num, dtype=np.float)
y = np.linspace(y_min, y_max, y_num, dtype=np.float)
t = np.arange(t_min, t_max+1, t_step, dtype=np.float)

tt, yy, xx = np.meshgrid(t, y, x, indexing="ij")


### Directories
current_dir = os.path.dirname(os.path.realpath(__file__))
pressure_dir = f"{current_dir}/pressure"

os.makedirs(pressure_dir, exist_ok=True)


### Function
def pressure(x, y, t, t0=10000., U=50., a=200000., p0=2000.):
    return (
        p0
        * (1. - np.exp(-t / t0))
        * np.exp(-(x**2. + (y - U * t)**2.) / a**2.)
    )


### Compute field
for case_number in range(num_cases):
    case = cases[case_number]
    U = U_list[case_number]
    a = a_list[case_number]
    p0 = p0_list[case_number]

    filename = f"{pressure_dir}/repr_{case:02.0f}"
    figurename = f"{pressure_dir}/test_repr_{case:02.0f}"

    print(
        f"\nComputing pressure field for case {case} (U={U}, a={a}, p0={p0})")

    # loop over all t, y and x
    # p = np.zeros((t.size, y.size, x.size))
    # for i in range(t.size):
    #     for j in range(y.size):
    #         for k in range(x.size):
    #             p[i,j,k] = pressure(x[k], y[j], t[i], t0, U, a, p0)

    p = pressure(xx, yy, tt, t0, U, a, p0)

    # remove zero-columns and zero-rows
    ix = np.where(~ np.all(np.isclose(p, 0), axis=(0, 1)))[0]
    iy = np.where(~ np.all(np.isclose(p, 0), axis=(0, 2)))[0]
    x = x[ix]
    y = y[iy]
    p = p[:, :, ix][:, iy, :]

    ### Write field
    print(f"Writing pressure field for case {case}")
    data = fp.convert_to_xarray(t, x, y, p)
    del t, x, y, p
    fp.write_pressure(data, filename)

    ### Write forcing file
    if case != 0:
        print(f"Overwriting forcing file for case {case}")
        with open(f"{current_dir}/forcing_repr_{case:02.0f}.ext", "w") as file:
            file.write("* Meteo forcing \n")
            file.write("QUANTITY = atmosphericpressure \n")
            file.write(f"FILENAME = pressure/repr_{case:02.0f}.amp \n")
            file.write("FILETYPE = 4 \n")
            file.write("METHOD   = 1 \n")
            file.write("OPERAND  = O \n")

    ### Visualise field
    print(f"Plotting pressure field for case {case}")
    fp.plot_pressure(data, filename=figurename, x_scales=[0, 500])


# plt.show()
