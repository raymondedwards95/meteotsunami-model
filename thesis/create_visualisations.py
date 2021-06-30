""" Script for making simple visualisations """

# fix for systems without display
import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys

import matplotlib.pyplot as plt
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.analysis as fa
import functions.animation as anim
import functions.regrid as fr
import functions.utilities as fu


### Parameters
parser = argparse.ArgumentParser(
    description="Process model outputs and create visualisations for a specific case"
)
parser.add_argument(
    "case",
    help="Case number",
    nargs="?",
    default="00",
    type=int
)
parser.add_argument(
    "--reprocess",
    "-r",
    help="Force processing of original data.",
    default=False,
    type=bool
)
parser.add_argument(
    "--show",
    "-s",
    help="Show figures after creation.",
    default=False,
    type=bool
)
parser.add_argument(
    "--animate",
    "-a",
    help="Create animation.",
    default=True,
    type=bool
)
args = parser.parse_args()

case = f"{args.case:02}"
reprocess_data = bool(args.reprocess)
show_figs = bool(args.show)
make_ani = bool(args.animate)


### Set paths
file_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = file_dir + "/output/"
output_dir = parent_dir + f"exp_{case}/"
figure_dir = file_dir + f"/figures/exp_{case}/"

filename_output = output_dir + "FlowFM_map.nc"
filename_processed = parent_dir + f"data_exp_{case}.nc"

filename_log = file_dir + "/last_runs.log"

print(f"Making figures for case {case}")
print(f"# Folder with data is \n# '{output_dir}'")
print(f"# File with raw data is \n# '{filename_output}'")
print(f"# File with processed data is \n# '{filename_processed}'")
print(f"# Folder that contains figures is \n# '{figure_dir}'")

# Check if data exists
if not os.path.exists(filename_output):
    print(f"File '{filename_output}' is not found. Exiting.")
    sys.exit(1)

# Convert data if it didn't happen yet
if not os.path.exists(filename_processed) or reprocess_data:
    print(f"File '{filename_processed}' is not found. Processing data.")
    fr.extract_data(filename_output, savename=filename_processed)
    with open(filename_log, "a+") as _file:
        _file.write(f"Processed data for case {case}")

os.makedirs(figure_dir, exist_ok=True)


### Get data
data = xr.open_dataset(filename_processed)


### Get optimal coordinates



### Figures


### Animation
if make_ani:
    anim.animation_alongshore(data, saveloc=figure_dir)
    anim.animation_contour(data, saveloc=figure_dir)
    anim.animation_contour_uv(data, saveloc=figure_dir)
else:
    print("\nNot creating animations")


### End
with open(filename_log, "a+") as _file:
    _file.write(f"Created visualisations for case {case}")
        
if show_figs: 
    plt.show()
