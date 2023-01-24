"""Script for creating animations from gridded model output """

import argparse
import os
import sys
import time

import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.animation as anim
import functions.utilities as fu
# fmt: on


# Parameters
t0 = time.perf_counter()

parser = argparse.ArgumentParser(
    description="Create animations for a specific exp case",
)
parser.add_argument(
    "case",
    help="Case number",
    nargs="?",
    default="00",
    type=int,
)
args = parser.parse_args()

case = f"{args.case:02}"

# Paths
current_dir = os.path.dirname(os.path.realpath(__file__))
anim_dir = f"{PATH_MAIN}/animations/exp_{case}"
data_dir = f"{current_dir}/output"

data_file = f"{data_dir}/data_exp_{case}.nc"

# Checks
print(f"Making animations for exp {case=}")
print(f"Save location: {anim_dir}")
print(f"Data file:     {data_file}")

if not os.path.exists(data_file):
    print(f"\nFile '{data_file}' is not found. Exiting.\n")
    sys.exit(1)

# Get data
data = xr.open_dataset(data_file)  # , chunks={"y": -1, "x": -1, "t": "auto"})
print(f"Grid size: {data.sizes}")
# print(f"Chunksize: {data.chunksizes}")

# Reduce data by skipping rows and columns with zeros
data_ref = np.max(np.abs(data["wl"].values))

ix = np.argwhere(
    ~np.all(np.abs(data["wl"].values) < 1e-2 * data_ref, axis=(0, 1))
).ravel()
iy = np.argwhere(
    ~np.all(np.abs(data["wl"].values) < 1e-2 * data_ref, axis=(0, 2))
).ravel()

x_min = np.floor(data["x"].values.min())
x_max = fu.relative_ceil(data["x"][np.max(ix)] / 5.0, s=1)

y_min = 0.0
y_max = fu.relative_ceil(data["y"][np.max(iy)], s=1)

data = data.sel(x=slice(x_min, x_max), y=slice(y_min, y_max))

# Downsample data
data = data.resample(t="600s").interpolate()
print(f"Downsampled grid size: {data.sizes}")
# print(f"Chunksize: {data.chunksizes}")

# Create animations
anim.animation_alongshore(data, savedir=anim_dir, fps=30.0)
anim.animation_contour(data, savedir=anim_dir, fps=30.0)
anim.animation_contour_uv(data, savedir=anim_dir, fps=30.0)

# Close data
data.close()

# End
t1 = time.perf_counter()
print(
    f"Finished creating animations for exp {case=} in {t1-t0:0.1f} seconds ({(t1-t0)/60.:0.1f} minutes)\n"
)
