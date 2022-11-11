"""Script for creating figures from gridded model output """

import argparse
import os
import sys
import time

import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.plotting as fpl
import functions.utilities as fu
# fmt: on


# Parameters
t0 = time.perf_counter()

parser = argparse.ArgumentParser(
    description="Create figures for a specific repr case",
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
figure_dir = f"{current_dir}/figures/repr_{case}"
data_dir = f"{current_dir}/output"

data_file = f"{data_dir}/data_repr_{case}.nc"

os.makedirs(figure_dir, exist_ok=True)

# Checks
print(f"Making figures for repr {case=}")
print(f"Save location: {figure_dir}")
print(f"Data file:     {data_file}")

if not os.path.exists(data_file):
    print(f"\nFile '{data_file}' is not found. Exiting.\n")
    sys.exit(1)

# Get data
data = xr.open_dataset(data_file)#, chunks={"y": -1, "x": -1, "t": "auto"})
print(f"Grid size: {data.sizes}")
# print(f"Chunksize: {data.chunksizes}")

# Parameters
t_list = np.linspace(
    fu.from_timestr(data["t"].min()),
    fu.from_timestr(data["t"].max()) - 7200.,
    5,)
y_list = np.linspace(
    0.,
    data["y"].max(),
    5,
)
y_min = 0.
x_max = np.round(data["x"].max() / 2. / 1000., -1)
x_ref = np.ceil(data["x"].min() * 2. / 1e3) * 1e3

# Contour
fpl.plot_contour() \
    .add_plots(
        dataset=data,
        variable_list=["wl"],
        t_list=t_list,
        y_min=y_min,
        x_max=x_max,
    ) \
    .save(figure_dir)

fpl.plot_contour() \
    .add_plots(
        dataset=data,
        variable_list=["wl", "p"],
        t_list=t_list,
        y_min=y_min,
        x_max=x_max,
) \
    .save(figure_dir)

fpl.plot_contour() \
    .add_plots(
        dataset=data,
        variable_list=["u", "v"],
        t_list=t_list,
        y_min=y_min,
        x_max=x_max,
) \
    .save(figure_dir)

# Alongshore
fpl.plot_alongshore(variable="wl") \
    .add_subplot(dataset=data,
                 x=x_ref, t=t_list,
    y_min=0,
    ) \
    .save(figure_dir)

# Crossshore
for t_single in t_list:
    fpl.plot_crossshore(variable="wl") \
        .plot_peaks(dataset=data, t=t_single, number=5,) \
            .save(figure_dir)

# Spectrum-1d
_spectrum_1d = fpl.plot_spectrum_1d(variable="wl", demean=True)
for y_single in y_list:
    _spectrum_1d.add_plot(dataset=data, x=x_ref, y=y_single,)
_spectrum_1d.save(figure_dir)
del _spectrum_1d

# Spectrum-2d
fpl.plot_spectrum_2d(variable="wl", demean=True) \
    .add_plot(dataset=data, x=x,) \
    .add_dispersion(n=3) \
    .save(figure_dir)

# Timeseries
_timeseries = fpl.plot_timeseries()
for y_single in y_list:
    _timeseries.add_plot(dataset=data, variable="wl", x=x_ref, y=y_single,)
_timeseries.save(figure_dir)
del _timeseries

# Close data
data.close()

# End
t1 = time.perf_counter()
print(
    f"Finished creating figures for repr {case=} in {t1-t0:0.1f} seconds ({(t1-t0)/60.:0.1f} minutes)\n"
)
