""" Script for making simple visualisations to compare with paper An et al. (2012) """
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.animation as anim
import functions.regrid as fr
import functions.utilities as fu
import functions.visualisation as fv


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
output_dir = parent_dir + f"repr_{case}/"
figure_dir = file_dir + f"/figures/repr_{case}/"

filename_output = output_dir + "FlowFM_map.nc"
filename_processed = parent_dir + f"data_repr_{case}.nc"

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
        _file.write(f"\tProcessed data for case {case}\n")

os.makedirs(figure_dir, exist_ok=True)


### Get data
print("Read data")
data = xr.open_dataset(filename_processed)
x = data["x"]
y = data["y"]
t = data["t"]

b = data["b"]
wl = data["wl"]
u = data["u"]
v = data["v"]
p = data["p"]

# offset x-coordinate for nicer plotting
x = x - x.min()


try:
    print(f"Creating figures in \n'{figure_dir}'")
    ### Figure bathymetry
    print("Figure: bathymetry")
    plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)

    plt.contourf(x/1000, y/1000, b)
    plt.colorbar()

    plt.title("Water depth [m]")
    plt.xlabel("$x$ [km]")
    plt.ylabel("$y$ [km]")

    plt.axhline(color="black", linewidth=1)
    plt.axvline(color="black", linewidth=1)

    # plt.axhline(7.56e3, linestyle="--", color="gray", linewidth=1)
    # plt.axvline(10, linestyle="--", color="gray", linewidth=1)
    # plt.axvline(100, linestyle="--", color="gray", linewidth=1)
    # plt.axvline(200, linestyle="--", color="gray", linewidth=1)
        
    plt.savefig(figure_dir + "bathymetry_contour", bbox_inches="tight", dpi=FIG_DPI)
except:
    print(f"Error in bottom-profile visualisation {case=}")


### Figure pressure
print("Figure: pressure")
plot_times = [4e4, 8e4, 12e4, 16e4]

try:
    fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
    fig.set_size_inches(FIGSIZE_SQUARE)
    fig.set_dpi(FIG_DPI)
    fig.suptitle("Pressure distribution [Pa]")
    fig.set_tight_layout(True)

    for i in range(4):
        _ax = ax[i//2, i%2]
        im = _ax.contourf(x/1000, y/1000, p.interp(t=fu.to_timestr(plot_times[i])), vmin=0, vmax=2000)
        _ax.set_xlim([0, 300])
        _ax.set_ylim([0, y.max()/1000])
        _ax.set_title(f"$t = {plot_times[i] : 0.0f}$s")
        
        if i//2: _ax.set_xlabel("$x$ [km]")
        if not i%2: _ax.set_ylabel("$y$ [km]")
        
        fig.colorbar(im, ax=_ax)
        
    #     _ax.axhline(7.56e3, linestyle="--", color="gray", linewidth=1)
    #     _ax.axvline(10, linestyle="--", color="gray", linewidth=1)
    #     _ax.axvline(100, linestyle="--", color="gray", linewidth=1)
    #     _ax.axvline(200, linestyle="--", color="gray", linewidth=1)
        
    fig.savefig(figure_dir + "pressure_contours", bbox_inches="tight", dpi=FIG_DPI)
except:
    print(f"Error in pressure visualisation {case=}")


### Figure sse contours
print("Figure: sse contours")
# figure 1
plot_times = [4e4, 8e4, 12e4, 16e4]

try:
    fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
    fig.set_size_inches(FIGSIZE_SQUARE)
    fig.set_dpi(FIG_DPI)
    fig.suptitle("Sea Surface Elevation [m]")
    for i in range(4):
        _ax = ax[i//2, i%2]
        im = _ax.contourf(x/1000, y/1000, wl.interp(t=fu.to_timestr(plot_times[i])), vmin=-0.8, vmax=0.8)
        _ax.set_xlim([0, 300])
        _ax.set_ylim([0, y.max()/1000])
        _ax.set_title(f"$t = {plot_times[i] : 0.0f}$s")
        
        if i//2: _ax.set_xlabel("$x$ [km]")
        if not i%2: _ax.set_ylabel("$y$ [km]")
        
        fig.colorbar(im, ax=_ax)
        
    #     _ax.axhline(7.56e3, linestyle="--", color="gray", linewidth=1)
    #     _ax.axvline(10, linestyle="--", color="gray", linewidth=1)
    #     _ax.axvline(100, linestyle="--", color="gray", linewidth=1)
    #     _ax.axvline(200, linestyle="--", color="gray", linewidth=1)
        
    fig.savefig(figure_dir + "sse_contours", bbox_inches="tight", dpi=FIG_DPI)
except:
    print(f"Error in contours visualisation {case=}")

### Figure alongshore profiles sse
# figure 2
print("Figure: sse alongshore profile")
plot_times = [4e4, 8e4, 12e4, 16e4]
plot_ylims = np.array([[0, 4], [1, 6], [2, 8], [3, 10]]) * 1e3

try:
    fig, ax = plt.subplots(2, 2, sharex=False, sharey=False)
    fig.set_size_inches(FIGSIZE_SQUARE)
    fig.set_dpi(FIG_DPI)
    fig.suptitle("Along-shore profile of Sea Surface Elevation [m]")
    for i in range(4):
        _ax = ax[i//2, i%2]
        im = _ax.plot(y/1000, wl.interp(t=fu.to_timestr(plot_times[i]), x=10e3))
        _ax.set_xlim(plot_ylims[i])
        _ax.set_ylim([-1, 1])
        _ax.set_title(f"$t = {plot_times[i] : 0.0f}$s")
        
        if i//2: _ax.set_xlabel("$y$ [km]")
        if not i%2: _ax.set_ylabel("$SSE$ [m]")
            
        _ax.axvline(7.56e3, linestyle="--", color="gray", linewidth=1)
        _ax.axhline(color="black", linewidth=1)
        _ax.axvline(color="black", linewidth=1)
        
    fig.savefig(figure_dir + "sse_along", bbox_inches="tight", dpi=FIG_DPI)
except:
    print(f"Error in alongshore-profile visualisation {case=}")


### Figure cross-shore profile sse
# figure 3
print("Figure: sse cross-shore profile")
y_slices = np.array([7.56]) * 1e6

try:
    plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)
    for i in range(len(y_slices)):
        plt.plot(x/1000, wl.interp(t=fu.to_timestr(1.6e5), y=y_slices[i]), label=f"$y={y_slices[i]/1000 : 0.0f}$ km")
    plt.legend()
    plt.axhline(color="black", linewidth=1)
    plt.axvline(color="black", linewidth=1)
    plt.title("Cross-shore profile of Sea Surface Elevation [m]")
    plt.xlabel("$x$ [km]")
    plt.ylabel("$SSE$ [m]")
    plt.xlim([0, 600])
    plt.ylim([0, 0.7])

    # plt.axvline(10, linestyle="--", color="gray", linewidth=1)
    # plt.axvline(100, linestyle="--", color="gray", linewidth=1)
    # plt.axvline(200, linestyle="--", color="gray", linewidth=1)
        
    plt.savefig(figure_dir + "sse_cross", bbox_inches="tight", dpi=FIG_DPI)
except:
    print(f"Error in cross-shore-profile visualisation {case=}")


### Figure alongshore profiles
# figure
print("Figure: sse more alongshore profiles")
x_slices = np.array([10, 100, 200]) * 1e3
linestyles = ["-", "--", "--"]

try:
    plt.figure(figsize=FIGSIZE_NORMAL, dpi=FIG_DPI)
    plt.title("Along-shore profile of Sea Surface Elevation [m]")
    for i in range(len(x_slices)):
        plt.plot(
            y/1000, 
            wl.interp(t=fu.to_timestr(1.6e5), x=x_slices[i]), 
            label=f"$x={x_slices[i]/1000 : 0.0f}$ km",
            linestyle=linestyles[i]
        )
    plt.legend()
    plt.axhline(color="black", linewidth=1)
    plt.axvline(color="black", linewidth=1)
    plt.xlim([0, y.max()/1000])
    plt.axvline(7.56e3, linestyle="--", color="gray", linewidth=1)
    plt.xlabel("$y$ [km]")
    plt.ylabel("$SSE$ [m]")
        
    plt.savefig(figure_dir + "sse_along_2", bbox_inches="tight", dpi=FIG_DPI)
except:
    print(f"Error in alongshore-profile 2 visualisation {case=}")


### Figures - Contour
t_moments = [4e4, 8e4, 12e4, 16e4]
try:
    fv.vis_contour(data, t=t_moments, saveloc=figure_dir, variable="wl", ylims=[0, 300])
    fv.vis_contour(data, t=t_moments, saveloc=figure_dir, variable="u", ylims=[0, 300])
    fv.vis_contour(data, t=t_moments, saveloc=figure_dir, variable="v", ylims=[0, 300])
except:
    print(f"Error in contour visualisation {case=}")


### Figures - Alongshore
try:
    fv.vis_alongshore(data, t=t_moments, saveloc=figure_dir)
    for t_ in t_moments:
        fv.vis_alongshore(data, t=t_, saveloc=figure_dir)
except:
    print(f"Error in alongshore visualisation {case=}")


### Figures - Timeseries
y_list_0 = np.array([4, 8]) * 1e6
y_list_1 = np.array([2, 4, 6, 8]) * 1e6
try:
    fv.vis_timeseries(data, y=y_list_0, saveloc=figure_dir)
    fv.vis_timeseries(data, y=y_list_1, saveloc=figure_dir)
except:
    print(f"Error in timeseries visualisation {case=}")


### Figures - Spectra 1d
try:
    for _y in np.union1d(y_list_0, y_list_1):
        fv.vis_spectrum_1d(data, x=1e4, y=_y, saveloc=figure_dir)
except:
    print(f"Error in spectrum-1d visualisation {case=}")
    
try:
    fv.vis_spectrum_1d(data, x=1e4, y=y_list_0, saveloc=figure_dir)
    fv.vis_spectrum_1d(data, x=1e4, y=y_list_1, saveloc=figure_dir)
except:
    print(f"Error in spectrum-1d visualisation {case=}")


### Figures - Spectra 2d
try:
    fv.vis_spectrum_2d(data, x=1e4, saveloc=figure_dir)
except:
    print(f"Error in spectrum-1d visualisation {case=}")


### Figures - Cross-shore
try:
    fv.vis_crossshore(data, y=y_list_0, t=t_moments, saveloc=figure_dir)
    fv.vis_crossshore(data, y=y_list_1, t=t_moments, saveloc=figure_dir)

    for t_ in t_moments:
        fv.vis_crossshore(data, y=y_list_0, t=t_, saveloc=figure_dir)
        fv.vis_crossshore(data, y=y_list_1, t=t_, saveloc=figure_dir)

    for y_ in np.union1d(y_list_0, y_list_1):
        fv.vis_crossshore(data, y=y_, t=t_moments, saveloc=figure_dir)

    for t_ in t_moments:
        for y_ in y_list_1:
            fv.vis_crossshore(data, y=y_, t=t_, saveloc=figure_dir)
except:
    print(f"Error in cross-shore visualisation {case=}")


### Animation
if make_ani:
    anim.animation_alongshore(data, saveloc=figure_dir)
    anim.animation_contour(data, saveloc=figure_dir)
    anim.animation_contour_uv(data, saveloc=figure_dir)
else:
    print("\nNot creating animations")


### End
with open(filename_log, "a+") as _file:
    _file.write(f"\tCreated visualisations for case {case}\n")
        
if show_figs: 
    plt.show()
