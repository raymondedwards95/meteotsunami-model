""" Creates a file with atmospheric pressure """

import os
import sys

import dask.array as da
import numpy as np

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.pressure as fp


### Directories
current_dir = os.path.dirname(os.path.realpath(__file__))
pressure_dir = f"{current_dir}/pressure"

os.makedirs(pressure_dir, exist_ok=True)


### Parameters
# pressure distribution
t0 = 10000.
U_list = np.array([5, 15, 25])
a_list = np.array([10000, 20000, 30000])
p0_list = np.array([2000])
x0_list = np.array([0, 50000])

# cross shore (meters)
x_min = 0
x_max = 1e6
x_step = 5e3

# along shore (meters)
y_min = -2e6
y_max = +3e6
y_step = x_step

# time (seconds)
t_min = 0
t_max = 50. * 3600.
t_step = 3600. / 6.


### Convert parameters
x0_list, a_list, p0_list, U_list = np.meshgrid(x0_list, a_list, p0_list, U_list, indexing="ij")
x0_list = np.ravel(x0_list)
a_list = np.ravel(a_list)
p0_list = np.ravel(p0_list)
U_list = np.ravel(U_list)

cases = np.arange(len(U_list))
num_cases = len(cases)


### Check parameters
assert len(U_list) == num_cases
assert len(a_list) == num_cases
assert len(p0_list) == num_cases
assert len(x0_list) == num_cases


### Save parameters to file
print("\nPressure fields for the following cases are computed (case, U, a, p0, x0)")
with open(f"{current_dir}/parameters_pressure.txt", "w") as file:
    for i in range(num_cases):
        line = f"{cases[i]:02.0f},{U_list[i]:0.0f},{a_list[i]:0.0f},{p0_list[i]:0.0f},{x0_list[i]:0.0f}"
        print(line.replace(",", "\t"))
        file.write(line + "\n")

    del line, i


### Grid
x_num = int((x_max - x_min) / x_step + 1)
y_num = int((y_max - y_min) / y_step + 1)

x = da.linspace(x_min, x_max, x_num, chunks=256)
y = da.linspace(y_min, y_max, y_num, chunks=1024)
t = da.arange(t_min, t_max+1, t_step, chunks=128)

t_num = t.size

tt, yy, xx = da.meshgrid(t, y, x, indexing="ij", sparse=True)

print("Grid parameters:")
print(f"{x.size=}\t\t{y.size=}\t\t{t.size=}")
print(f"{x.chunksize=}\t{y.chunksize=}\t{t.chunksize=}")


### Function
def pressure(x, y, t, t0=10000., U=50., a=200000., p0=2000., x0=0.):
    return (
        p0
        * (1. - da.exp(-t / t0))
        * da.exp(-((x - x0)**2. + (y - U * t)**2.) / a**2.)
    )


### Compute fields
for case_number in range(num_cases):
    ## Set parameters
    case = cases[case_number]
    U = U_list[case_number]
    a = a_list[case_number]
    p0 = p0_list[case_number]
    x0 = x0_list[case_number]

    ## Set paths
    filename = f"{pressure_dir}/exp_{case:02.0f}"
    figurename = f"{pressure_dir}/fig_exp_{case:02.0f}"

    ## Compute pressure
    print(f"\nComputing pressure field for {case=} ({U=}, {a=}, {p0=}, {x0=})")
    p = pressure(xx, yy, tt, t0, U, a, p0, x0).astype(np.float32)

    ## Remove zero-columns and zero-rows
    ix = da.where(~ da.all(da.isclose(p, 0), axis=(0, 1)))[0]
    iy = da.where(~ da.all(da.isclose(p, 0), axis=(0, 2)))[0]
    _x = x[ix]
    _y = y[iy]
    p = p[:, :, ix][:, iy, :]

    ## Write field
    data = fp.convert_to_xarray(t, _x, _y, p.compute())
    print(f"Process pressure field for {case=}")
    del _x, _y, p
    print(f"Writing pressure field for {case=}")
    fp.write_pressure(data, filename)

    ## Write forcing file
    print(f"Overwriting forcing file for {case=}")
    with open(f"{current_dir}/pressure/forcing_exp_{case:02.0f}.ext", "w") as file:
        file.write("* Meteo forcing \n")
        file.write("QUANTITY = atmosphericpressure \n")
        file.write(f"FILENAME = exp_{case:02.0f}.amp \n")
        file.write("FILETYPE = 4 \n")
        file.write("METHOD   = 1 \n")
        file.write("OPERAND  = O \n")

    ## Visualise field
    print(f"Plotting pressure field for {case=}")
    fp.plot_pressure(data, filename=figurename)


# plt.show()
print("Finished creating pressure-files")
