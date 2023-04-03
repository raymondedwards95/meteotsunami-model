""" Creates a file for D3D-FM-FLOW with a bathymetry as in paper An et al. (2012) """

import os
import sys
import time

import numpy as np

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.bathymetry as fb
# fmt: on

t0 = time.perf_counter()
print(f"\nStart creating bathymetry-files for repr")


# Parameters
cases = [0, 36, 41, 42, 43]
slopes = [1.0 / 400.0, 1.0 / 800.0, 0.0, 0.0, 0.0]
depths = [0, 0, 250, 100, 500]


# Grid
x, dx = np.linspace(0, 1e6, 7, retstep=True)
y, dy = np.linspace(-1e7, 1e7, 7, retstep=True)
xx, yy = np.meshgrid(x, y)


# Function
def simple_slope(x, slope=1.0 / 400.0, ref_depth=0.0):
    """Computes bed height as `b = - ref_depth - slope * x`"""
    return -1.0 * (ref_depth + slope * x)


# Paths
script_dir = os.path.dirname(os.path.realpath(__file__))
bathymetry_dir = f"{script_dir}/bathymetry"
os.makedirs(bathymetry_dir, exist_ok=True)


# Save parameters to file
print("\nBathymetries for the following cases are computed: \ncase \tslope \tdepth")
with open(f"{bathymetry_dir}/parameters_bathymetry.txt", "w") as file:
    file.write(f"Sloped bottom parameters\n")
    file.write(f"case,slope,depth\n")
    for i in range(len(cases)):
        line = f"{cases[i]:02.0f},{slopes[i]:0.2e},{depths[i]:0.0f}"
        line = line.replace("0.00e+00", "0")
        print(line.replace(",", "\t"))
        file.write(f"{line}\n")

    file.write(f"\nGrid parameters\n")
    file.write(f"coordinate,min,max,num,step\n")
    file.write(f"x,{x.min()},{x.max()},{x.size},{dx:0.0f}\n")
    file.write(f"y,{y.min()},{y.max()},{y.size},{dy:0.0f}\n")

    del line, i


# Computations
for i in range(len(cases)):
    # Start
    ta = time.perf_counter()

    # Take parameters
    case = cases[i]
    slope = slopes[i]
    depth_0 = depths[i]
    print(f"\n### Creating bathymetry for {case=:02.0f}, with {slope=} and {depth_0=}")

    # Compute bathymetry
    zz = simple_slope(xx, slope, depth_0)

    # Write to file
    data = fb.convert_to_xarray(
        x,
        y,
        zz,
        savename=f"{bathymetry_dir}/repr_{case:02.0f}.nc",
    )
    fb.write_bathymetry(
        data,
        f"{bathymetry_dir}/repr_{case:02.0f}.xyb",
    )

    # Visualize
    fb.plot_bathymetry(
        data,
        f"{bathymetry_dir}/repr_{case:02.0f}",
    )
    fb.plot_bathymetry(
        data,
        f"{bathymetry_dir}/repr_{case:02.0f}",
        half_width=True,
    )

    # End
    tb = time.perf_counter()
    print(f"### Finished creating bathymetry for {case=:02.0f} in {tb-ta:0.1f} seconds")


# End
t1 = time.perf_counter()
print(f"\nFinished creating bathymetry-files for repr in {t1-t0:0.1f} seconds")
