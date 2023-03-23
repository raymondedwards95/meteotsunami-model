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
figure_dir = f"{current_dir}/results/repr_{case}"
data_dir = f"{current_dir}/output"

data_file = f"{data_dir}/data_repr_{case}.nc"

# Checks
print(f"Making figures for repr {case=}")
print(f"Save location: {figure_dir}")
print(f"Data file:     {data_file}")

if not os.path.exists(data_file):
    print(f"\nFile '{data_file}' is not found. Exiting.\n")
    sys.exit(1)

# Get data
data = xr.open_dataset(data_file, chunks={"y": -1, "x": -1, "t": "auto"})
print(f"Grid size: {data.sizes}")
print(f"Chunksize: {data.chunksizes}")

# Parameters
t_list = np.linspace(
    fu.from_timestr(data["t"].values.min()),
    fu.from_timestr(data["t"].values.max()) - 7200.0,
    5,
)
y_list = np.linspace(
    0.0,
    data["y"].max(),
    5,
)
y_min = 0.0
x_max = np.round(data["x"].values.max() / 2.0 / 1000.0, -1)
x_ref = np.ceil(data["x"].values.min() * 2.0 / 1e3) * 1e3

print(f"{x_max=}")
print(f"{x_ref=}")

# Find limits
data_ref = np.max(np.abs(data["wl"].values))

ix = np.argwhere(
    ~np.all(np.abs(data["wl"].values) < 1e-2 * data_ref, axis=(0, 1))
).ravel()
iy = np.argwhere(
    ~np.all(np.abs(data["wl"].values) < 1e-2 * data_ref, axis=(0, 2))
).ravel()

x_min = np.floor(data["x"].values.min())
x_max = fu.relative_ceil(data["x"][np.max(ix)] / 2.0, s=2)

y_min = 0.0
y_max = fu.relative_ceil(data["y"][np.max(iy)], s=2)

# print(f"{x_ref=}")
# print(f"{x_min=}")
# print(f"{x_max=}")
# print(f"{y_min=}")
# print(f"{y_max=}")

# Contour
_contour_a = fpl.plot_contour(scale="Mm")
_contour_a.add_plots(
    dataset=data,
    variable_list=["p"],
    t_list=t_list,
    y_min=y_min,
    y_max=y_max,
    x_max=x_max,
)
_contour_a.save(figure_dir)

_contour_b = fpl.plot_contour(scale="Mm")
_contour_b.add_plots(
    dataset=data,
    variable_list=["wl"],
    t_list=t_list,
    y_min=y_min,
    y_max=y_max,
    x_max=x_max,
)
_contour_b.save(figure_dir)

_contour_c = fpl.plot_contour(scale="Mm")
_contour_c.add_plots(
    dataset=data,
    variable_list=["wl", "p"],
    t_list=t_list,
    y_min=y_min,
    y_max=y_max,
    x_max=x_max,
)
_contour_c.save(figure_dir)

_contour_d = fpl.plot_contour(scale="Mm")
_contour_d.add_plots(
    dataset=data,
    variable_list=["u", "v"],
    t_list=t_list,
    y_min=y_min,
    y_max=y_max,
    x_max=x_max,
)
_contour_d.save(figure_dir)

# Alongshore
_along_shore = fpl.plot_alongshore(variable="wl", scale="Mm")
_along_shore.add_subplot(
    dataset=data,
    x=x_ref,
    t=t_list,
    y_min=0,
    y_max=y_max,
)
_along_shore.save(figure_dir)

# Crossshore
for t_single in t_list:
    _cross_shore = fpl.plot_crossshore(variable="wl", scale="Mm")
    _cross_shore.plot_peaks(
        dataset=data,
        t=t_single,
        number=5,
        x_max=x_max,
    )
    _cross_shore.save(figure_dir)

# Spectrum-1d
_spectrum_1d = fpl.plot_spectrum_1d(variable="wl", demean=True)
for y_single in y_list[1:-1]:
    _spectrum_1d.add_plot(
        dataset=data,
        x=x_ref,
        y=y_single,
        label=f"\\( y = \\SI{{{y_single / 1e3:0.1f}}}{{\\kilo\\meter}} \\)",
    )
_spectrum_1d.save(figure_dir)

# Spectrum-2d
_spectrum_2d = fpl.plot_spectrum_2d(variable="wl", demean=True, scale="Mm")
_spectrum_2d.add_plot(
    dataset=data,
    x=x_ref,
)
_spectrum_2d.add_dispersion(n=3)
_spectrum_2d.save(figure_dir)

# Timeseries
_timeseries = fpl.plot_timeseries()
for y_single in y_list:
    _timeseries.add_plot(
        dataset=data,
        variable="wl",
        x=x_ref,
        y=y_single,
    )
_timeseries.save(figure_dir)

# Close data
data.close()

# End
t1 = time.perf_counter()
print(
    f"Finished creating figures for repr {case=} in {t1-t0:0.1f} seconds ({(t1-t0)/60.:0.1f} minutes)\n"
)
