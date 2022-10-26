""" Creates files to describe atmospheric pressure for D3D-FM-FLOW for the experiments """

import os
import sys
import time

import dask.array as da
import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp
# fmt: on


# Start
t0 = time.perf_counter()
print(f"\nStart creating bathymetry-files for exp")


# Parameters
def pressure(x, y, t, t0=10000.0, U=50.0, a=200000.0, p0=2000.0, x0=0.0):
    """Pressure disturbance distribution used for experiments

    Input:
        x:  array of x-coordinates
        y:  array of y-coordinates
        t:  array of time-coordinates

    Options:
        t0: growth-timescale factor
        U:  travel velocity of pressure disturbance
        a:  size of pressure disturbance
        p0: magnitude of pressure disturbance
        x0: x-coordinate of the center of the pressure disturbance

    Output:
        p:  pressure
    """
    return (
        p0
        * (1.0 - da.exp(-t / t0))
        * da.exp(-((x - x0) ** 2.0 + (y - U * t) ** 2.0) / a**2.0)
    )


# pressure distribution
t0_value = 10000
U_list = np.array([5, 15, 25], dtype=np.float32)
a_list = np.array([10000, 20000, 30000], dtype=np.float32)
p0_list = np.array([2000], dtype=np.float32)
x0_list = np.array([0, 50000], dtype=np.float32)

# cross shore (meters)
x_min = 0
x_max = 1e6
x_step = 2.5e3

# along shore (meters)
y_min = -2e6
y_max = +3e6
y_step = x_step

# time (seconds)
t_min = 0
t_max = 50.0 * 3600.0
t_step = 3600.0 / 10.0


# Directories
current_dir = os.path.dirname(os.path.realpath(__file__))
pressure_dir = f"{current_dir}/pressure"

os.makedirs(pressure_dir, exist_ok=True)


# Convert parameters
x0_array, a_array, p0_array, U_array = np.meshgrid(
    x0_list,
    a_list,
    p0_list,
    U_list,
    indexing="ij",
)
x0_array = np.ravel(x0_array)
a_array = np.ravel(a_array)
p0_array = np.ravel(p0_array)
U_array = np.ravel(U_array)

cases = np.arange(len(U_array))
num_cases = np.size(cases)


# Check parameters
assert np.size(U_array) == num_cases
assert np.size(a_array) == num_cases
assert np.size(p0_array) == num_cases
assert np.size(x0_array) == num_cases


# Save parameters to file
print(
    "\nPressure fields for the following cases are computed: \ncase \tU \ta \tp0 \tx0"
)
with open(f"{pressure_dir}/parameters_pressure.txt", "w") as file:
    file.write(f"case,U,a,p0,x0\n")
    for i in range(num_cases):
        line = f"{cases[i]:02.0f},{U_array[i]:0.0f},{a_array[i]:0.0f},{p0_array[i]:0.0f},{x0_array[i]:0.0f}"
        print(line.replace(",", "\t"))
        file.write(f"{line}\n")

    del line, i


# Grid
x_num = int((x_max - x_min) / x_step + 1)
y_num = int((y_max - y_min) / y_step + 1)

x = da.linspace(x_min, x_max, x_num, chunks="auto", dtype=np.float32)
y = da.linspace(y_min, y_max, y_num, chunks="auto", dtype=np.float32)
t = da.arange(t_min, t_max + 1, t_step, chunks=20, dtype=np.float32)

t_num = t.size

tt, yy, xx = da.meshgrid(t, y, x, indexing="ij", sparse=True)

print("Grid parameters:")
print(f"{x.size=}\t\t{y.size=}\t\t{t.size=}")
print(f"{x.chunksize=}\t{y.chunksize=}\t{t.chunksize=}")


# Compute fields
for case_number in range(num_cases):
    # Set parameters
    case = cases[case_number]
    U = U_array[case_number]
    a = a_array[case_number]
    p0 = p0_array[case_number]
    x0 = x0_array[case_number]

    # Set paths
    filename = f"{pressure_dir}/exp_{case:02.0f}"
    figurename = f"{pressure_dir}/fig_exp_{case:02.0f}"

    # Compute pressure
    print(f"\nComputing pressure field for {case=:02.0f} ({U=}, {a=}, {p0=}, {x0=})")
    p = pressure(xx, yy, tt, t0_value, U, a, p0, x0).astype(np.float32)

    # Process pressure
    print(f"Process pressure field for {case=:02.0f}")
    fp.convert_to_xarray(t, x, y, p, savename=filename, close=True)
    del p

    # Re-read data (lazy)
    print(f"Read data for {case=:02.0f}")
    data = xr.open_dataarray(f"{filename}.nc", chunks="auto")

    # Write field
    print(f"Writing pressure field for {case=:02.0f}")
    fp.write_pressure(data, filename)

    # Write forcing file
    print(f"Overwriting forcing file for {case=:02.0f}")
    with open(f"{pressure_dir}/forcing_exp_{case:02.0f}.ext", "w") as file:
        file.write("* Meteo forcing \n")
        file.write("QUANTITY = atmosphericpressure \n")
        file.write(f"FILENAME = exp_{case:02.0f}.amp \n")
        file.write("FILETYPE = 4 \n")
        file.write("METHOD   = 1 \n")
        file.write("OPERAND  = O \n")

    # Visualise field
    print(f"Plotting pressure field for {case=:02.0f}")
    fp.plot_pressure(data, filename=figurename)


# plt.show()
t1 = time.perf_counter()
print(f"\nFinished creating pressure-files for exp in {t1-t0:0.1f} seconds")
