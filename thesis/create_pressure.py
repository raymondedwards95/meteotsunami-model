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
cases = [0, 1]
num_cases = len(cases)

# pressure distribution
t0 = 10000.
U_list = [50.] * 2
a_list = [20000.] * 2
p0_list = [2000] * 2 
x0_list = [0., 5e4]

# cross shore (meters)
x_min = 0.  
x_max = 1e5
x_step = 10e3  

# along shore (meters)
y_min = -1e6
y_max = 1e6
y_step = x_step  

# time (seconds)
t_min = 0
t_max = 70 * 3600.  
t_step = 3600. / 4.  

# check parameters
assert len(U_list) == num_cases
assert len(a_list) == num_cases
assert len(p0_list) == num_cases
assert len(x0_list) == num_cases

print("\nPressure fields for the following cases are computed (case, U, a, p0, x0)")
for i in range(num_cases):
    print(cases[i], U_list[i], a_list[i], p0_list[i], x0_list[i])


### Grid
x_num = int((x_max - x_min) / x_step + 1)
y_num = int((y_max - y_min) / y_step + 1)

x = np.linspace(x_min, x_max, x_num)
y = np.linspace(y_min, y_max, y_num)
t = np.arange(t_min, t_max+1, t_step)

tt, yy, xx = np.meshgrid(t, y, x, indexing="ij")


### Directories
current_dir = os.path.dirname(os.path.realpath(__file__))
pressure_dir = f"{current_dir}/pressure"

os.makedirs(pressure_dir, exist_ok=True)


### Function
def pressure(x, y, t, t0=10000., U=50., a=200000., p0=2000., x0=0.):
    return (
        p0
        * (1. - np.exp(-t / t0))
        * np.exp(-((x - x0)**2. + (y - U * t)**2.) / a**2.)
    )


### Compute field
for case_number in range(num_cases):
    case = cases[case_number]
    U = U_list[case_number]
    a = a_list[case_number]
    p0 = p0_list[case_number]
    x0 = x0_list[case_number]

    filename = f"{pressure_dir}/exp_{case:02.0f}"
    figurename = f"{pressure_dir}/fig_exp_{case:02.0f}"

    print(
        f"\nComputing pressure field for case {case} (U={U}, a={a}, p0={p0})")

    # loop over all t, y and x
    # p = np.zeros((t.size, y.size, x.size))
    # for i in range(t.size):
    #     for j in range(y.size):
    #         for k in range(x.size):
    #             p[i,j,k] = pressure(x[k], y[j], t[i], t0, U, a, p0)

    p = pressure(xx, yy, tt, t0, U, a, p0, x0)

    # remove zero-columns and zero-rows
    ix = np.where(~ np.all(np.isclose(p, 0), axis=(0, 1)))[0]
    iy = np.where(~ np.all(np.isclose(p, 0), axis=(0, 2)))[0]
    _x = x[ix]
    _y = y[iy]
    p = p[:, :, ix][:, iy, :]

    ### Write field
    print(f"Writing pressure field for case {case}")
    data = fp.convert_to_xarray(t, _x, _y, p)
    del _x, _y, p
    fp.write_pressure(data, filename)

    ### Write forcing file
    if case != 0:
        print(f"Overwriting forcing file for case {case}")
        with open(f"{current_dir}/forcing_exp_{case:02.0f}.ext", "w") as file:
            file.write("* Meteo forcing \n")
            file.write("QUANTITY = atmosphericpressure \n")
            file.write(f"FILENAME = pressure/exp_{case:02.0f}.amp \n")
            file.write("FILETYPE = 4 \n")
            file.write("METHOD   = 1 \n")
            file.write("OPERAND  = O \n")

    ### Visualise field
    print(f"Plotting pressure field for case {case}")
    fp.plot_pressure(data, filename=figurename)


# plt.show()
print("Finished creating pressure-files")
