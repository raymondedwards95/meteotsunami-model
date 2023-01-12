""" Creates files for D3D-FM-FLOW with the pressure-field as in the paper An et al. (2012) """

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
print(f"\nStart creating bathymetry-files for repr")


# Parameters
# generic
cases = [0, 10, 11, 12, 15, 16, 17, 31, 32, 33, 41]
num_cases = len(cases)

# pressure distribution
# defaults:
# T0: 10000 s
# U: 50 m/s
# a: 200 km
# p0: 2000 Pa
T0 = 10000.0
U = 50.0
a = 200000.0
p0 = 2000.0

# cross shore (meters)
# x_min: 0 km
# x_max: 1000 km
# x_steps: 10 km, i.e. same as model-grid
x_min = 0.0
x_max = 1e6
x_steps = (
    np.array(
        [10, 5, 20, 40, 10, 10, 10, 10, 10, 10, 10],
        dtype=np.float32,
    )
    * 1e3
)

# along shore (meters)
# y_min: -10000 km
# y_max: 10000 km
# y_steps: 10 km, i.e. same as model-grid
y_min = -1e7
y_max = 1e7
y_steps = x_steps

# time (seconds)
# defaults:
# t_min: 0 hours
# t_max: 70 hours
# t_step: 30 minutes
t_min = 0
t_max = 55 * 3600.0
t_steps = (
    np.array(
        [30, 30, 30, 30, 60, 20, 15, 30, 30, 30, 30],
        dtype=np.float32,
    )
    * 60.0
)

# other
# x0 is the x-coordinate of the center
x0_vals = np.full(num_cases, 0, dtype=np.float32)
x0_vals[7] = 1e5
x0_vals[8] = 2e5
x0_vals[9] = 3e5
x0_vals[10] = 5e5


# Directories
current_dir = os.path.dirname(os.path.realpath(__file__))
pressure_dir = f"{current_dir}/pressure"

os.makedirs(pressure_dir, exist_ok=True)


# Check parameters
assert len(x_steps) == num_cases
assert len(y_steps) == num_cases
assert len(t_steps) == num_cases
assert len(x0_vals) == num_cases

print(
    "\nPressure fields for the following cases are computed: \ncase \tx_step \ty_step \tt_step \tx0"
)
with open(f"{pressure_dir}/parameters_pressure.txt", "w") as file:
    file.write(f"case,x_step,y_step,t_step,x0\n")
    for i in range(num_cases):
        line = f"{cases[i]:02.0f},{x_steps[i]:0.0f},{y_steps[i]:0.0f},{t_steps[i]:0.0f},{x0_vals[i]:0.0f}"
        print(line.replace(",", "\t"))
        file.write(f"{line}\n")

    del line, i


# Function
def pressure(x, y, t, t0=T0, U=U, a=a, p0=p0, x0=0.0):
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


# Compute field
for case_number in range(num_cases):
    # Start
    ta = time.perf_counter()

    # Set parameters
    case = cases[case_number]
    x_step = x_steps[case_number]
    y_step = y_steps[case_number]
    t_step = t_steps[case_number]
    x0 = x0_vals[case_number]

    print(
        f"\nComputing pressure field for case {case:02.0f} ({x_step:0.1f}, {y_step:0.1f}, {t_step:0.1f}, {x0:0.1f})"
    )

    # Set paths
    filename = f"{pressure_dir}/repr_{case:02.0f}"
    figurename = f"{pressure_dir}/repr_{case:02.0f}"

    # Create grid
    x_num = int((x_max - x_min) / x_step + 1)
    y_num = int((y_max - y_min) / y_step + 1)

    x, dx = da.linspace(x_min, x_max, x_num, chunks=-1, dtype=np.float32, retstep=True)
    y, dy = da.linspace(y_min, y_max, y_num, chunks=-1, dtype=np.float32, retstep=True)
    t = da.arange(t_min, t_max + 1, t_step, chunks=-1, dtype=np.float32)
    dt = t_step

    tt, yy, xx = da.meshgrid(t, y, x, indexing="ij")
    tt = tt.rechunk("auto")
    yy = yy.rechunk(tt.chunksize)
    xx = xx.rechunk(tt.chunksize)

    print("Grid parameters:")
    print(f"{t.size=}\t\t{y.size=}\t\t{x.size=}")
    print(f"{dt=:0.1f}\t\t{dy=:0.1f}\t\t{dx=:0.1f}")
    print(f"{tt.chunksize=}")
    print(f"{xx.chunksize=}")
    print(f"{yy.chunksize=}")

    # Compute pressure
    print(f"Computing pressure for case {case:02.0f}")
    p = pressure(xx, yy, tt, T0, U, a, p0, x0).astype(np.float32)

    # Convert data
    print(f"Converting data for case {case:02.0f}")
    fp.convert_to_xarray(t, x, y, p, savename=filename, close=True)
    del t, x, y, p

    # Re-read data (lazy)
    print(f"Read data for case {case:02.0f}")
    data = xr.open_dataarray(f"{filename}.nc", chunks="auto")

    # Write pressure data
    print(f"Writing pressure field for case {case:02.0f}")
    fp.write_pressure(data, filename)

    # Write forcing file
    if case != 0:
        print(f"Overwriting forcing file for case {case:02.0f}")
        with open(f"{current_dir}/pressure/forcing_repr_{case:02.0f}.ext", "w") as file:
            file.write("* Meteo forcing \n")
            file.write("QUANTITY = atmosphericpressure \n")
            file.write(f"FILENAME = repr_{case:02.0f}.amp \n")
            file.write("FILETYPE = 4 \n")
            file.write("METHOD   = 1 \n")
            file.write("OPERAND  = O \n")

    # Visualise field
    print(f"Plotting pressure field for case {case:02.0f}")
    fp.plot_pressure(data, filename=figurename, scale="Mm")

    # End
    tb = time.perf_counter()
    print(f"Finished creating pressure-field for {case=:02.0f} in {tb-ta:0.1f} seconds")


# End
t1 = time.perf_counter()
print(f"\nFinished creating pressure-files for repr in {t1-t0:0.1f} seconds")
