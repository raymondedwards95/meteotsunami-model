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
# fmt: on


# Parameters
t0 = time.perf_counter()

parser = argparse.ArgumentParser(
    description="Create animations for a specific repr case",
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
anim_dir = f"{PATH_MAIN}/animations/repr_{case}"
data_dir = f"{current_dir}/output"

data_file = f"{data_dir}/data_repr_{case}.nc"

# Checks
print(f"Making animations for repr {case=}")
print(f"Save location: {anim_dir}")
print(f"Data file:     {data_file}")

if not os.path.exists(data_file):
    print(f"\nFile '{data_file}' is not found. Exiting.\n")
    sys.exit(1)

# Get data
data = xr.open_dataset(data_file)  # , chunks={"y": -1, "x": -1, "t": "auto"})
print(f"Initial grid size:     {data.sizes}")
# print(f"Chunksize: {data.chunksizes}")

# Downsample data
data = data.resample(t="1800s").interpolate()
print(f"Downsampled grid size: {data.sizes}")
# print(f"Chunksize: {data.chunksizes}")

# Create animations
anim.animation_alongshore(data, savedir=anim_dir, fps=10.0)
anim.animation_contour(data, savedir=anim_dir, fps=10.0)
anim.animation_contour_uv(data, savedir=anim_dir, fps=10.0)

# Close data
data.close()

# End
t1 = time.perf_counter()
print(
    f"Finished creating animations for repr {case=} in {t1-t0:0.1f} seconds ({(t1-t0)/60.:0.1f} minutes)\n"
)
