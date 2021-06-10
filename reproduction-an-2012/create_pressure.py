""" Creates a file for D3D-FM-FLOW with the pressure-field as in the paper An et al. (2012) """

# fix for systems without display
import matplotlib
matplotlib.use("Agg")

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp


### Parameters
# generic
cases = [0, 10, 11, 12, 15, 16, 17]
num_cases = len(cases)

# pressure distribution
t0 = 10000.  # default: 10000 s
U = 50.  # default: 50 m/s
a = 200000.  # default: 200 km
p0 = 2000.  # default: 2000 Pa

# cross shore (meters)
x_min = 0.  # default: 0 km
x_max = 1e6  # default: 1000 km
x_steps = [1e4, 2e4, 4e4, 0.5e4, 1e4, 1e4, 1e4]  # default: 10 km

# along shore (meters)
y_min = -1e7  # default: -10000 km
y_max = 1e7  # default: 10000 km
y_steps = x_steps  # default: 10 km

# time (seconds)
t_min = 0  # default: 0
t_max = 70 * 3600.  # default: 70 hours
t_steps = list(np.array([1., 1., 1., 1., 2., 1./2., 1./4.]) * 3600.)  # default: 1 hour

# check parameters
assert len(x_steps) == num_cases
assert len(y_steps) == num_cases
assert len(t_steps) == num_cases

print("\nPressure fields for the following cases are computed (case, x_step, y_step, t_step)")
for i in range(num_cases):
    print(cases[i], x_steps[i], y_steps[i], t_steps[i])


### Directories
current_dir = os.path.dirname(os.path.realpath(__file__))
pressure_dir = f"{current_dir}/pressure"

os.makedirs(pressure_dir, exist_ok=True)


### Function
def pressure(x, y, t, t0=t0, U=U, a=a, p0=p0):
    return (
        p0
        * (1. - np.exp(-t / t0))
        * np.exp(-(x**2. + (y - U * t)**2. ) / a**2.)
    )


### Compute field
for case_number in range(num_cases):
    case = cases[case_number]
    x_step = x_steps[case_number]
    y_step = y_steps[case_number]
    t_step = t_steps[case_number]

    filename = f"{pressure_dir}/repr_{case:02.0f}"
    figurename = f"{pressure_dir}/fig_repr_{case:02.0f}"

    x_num = int((x_max - x_min) / x_step + 1)
    y_num = int((y_max - y_min) / y_step + 1)

    print(f"\nComputing pressure field for case {case} ({x_step}, {y_step}, {t_step})")
    x = np.linspace(x_min, x_max, x_num)
    y = np.linspace(y_min, y_max, y_num)
    t = np.arange(t_min, t_max+1, t_step)

    # loop over all t, y and x
    # p = np.zeros((t.size, y.size, x.size))
    # for i in range(t.size):
    #     for j in range(y.size):
    #         for k in range(x.size):
    #             p[i,j,k] = pressure(x[k], y[j], t[i], t0, U, a, p0)

    tt, yy, xx = np.meshgrid(t, y, x, indexing="ij")
    p = pressure(xx, yy, tt, t0, U, a, p0)

    # remove zero-columns and zero-rows
    ix = np.where(~ np.all(np.isclose(p, 0), axis=(0,1)))[0]
    iy = np.where(~ np.all(np.isclose(p, 0), axis=(0,2)))[0]
    x = x[ix]
    y = y[iy]
    p = p[:,:,ix][:,iy,:]


    ### Write field
    print(f"Writing pressure field for case {case}")
    data = fp.convert_to_xarray(t, x, y, p)
    del t, x, y, p
    fp.write_pressure(data, filename)


    ### Write forcing file
    if case != 0:
        print(f"Overwriting forcing file for case {case}")
        with open(f"{current_dir}/pressure/forcing_repr_{case:02.0f}.ext", "w") as file:
            file.write("* Meteo forcing \n")
            file.write("QUANTITY = atmosphericpressure \n")
            file.write(f"FILENAME = repr_{case:02.0f}.amp \n")
            file.write("FILETYPE = 4 \n")
            file.write("METHOD   = 1 \n")
            file.write("OPERAND  = O \n")


    ### Visualise field
    print(f"Plotting pressure field for case {case}")
    fp.plot_pressure(data, filename=figurename, x_scales=[0, 500])


# plt.show()
print("Finished creating pressure-files")
