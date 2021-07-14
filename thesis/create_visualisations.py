""" Script for making simple visualisations """
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.analysis as fa
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
        _file.write(f"Processed data for case {case}\n")

os.makedirs(figure_dir, exist_ok=True)


### Get data
data = xr.open_dataset(filename_processed)

t_moments = np.array([8, 16, 24, 32, 40]) * 3600.
x_moment = 1e4
x_offset = 5e4
y_moment = 500000.

### Compute stuff
with open(f"{figure_dir}/computed_parameters.txt", "w") as file:
    try:
        file.write(f"\n\nWave Period at y={y_moment} and x={x_moment}:\n")
        waveperiods = fa.compute_wave_periods(data, y_moment, x=x_moment)
        for waveperiod in waveperiods:
            file.write(f"{waveperiod:0.1f}\t")
        file.write(f"\nMean Wave Period:\n{np.nanmean(waveperiods)}\n")
    except:
        print(f"Error in computing waveperiod {case=}")


    for j in range(t_moments.size):
        t_moment = t_moments[j]

        try:
            file.write(f"\n\nWave Lengths at t={t_moment} and x={x_moment}:\n")
            wavelengths = fa.compute_wave_lengths(data, t_moment, x=x_moment)
            for wavelength in wavelengths:
                file.write(f"{wavelength:0.1f}\t")
            file.write(f"\nMean Wave Length:\n{np.nanmean(wavelengths)}")
            file.write(f"\nMean Wave Speed = Length / Period\n{np.nanmean(wavelengths) / np.nanmean(waveperiods)}\n")
        except:
            print(f"Error in computing wavelength {case=}")


### Figures
for j in range(t_moments.size):
    t_moment = t_moments[j]

    try:
        ## Find best parameters
        y_idx_max = fu.find_peaks_const_t(data, t_moment, x=x_moment, crests=True)
        y_idx_min = fu.find_peaks_const_t(data, t_moment, x=x_moment, crests=False)

        ## Cross-shore
        for i in range(np.min([y_idx_max.size, y_idx_min.size, 5])):
            fv.vis_crossshore(data, y=data["y"][y_idx_max[i]].values, t=t_moment, saveloc=figure_dir)
            fv.vis_crossshore(data, y=data["y"][y_idx_min[i]].values, t=t_moment, saveloc=figure_dir)
    except:
        print(f"Error in cross-shore visualisation {case=}")

    ## Along-shore
    try:
        fv.vis_alongshore(data, t=t_moment, x=x_moment, saveloc=figure_dir)
        fv.vis_alongshore(data, t=t_moment, x=x_moment+x_offset, saveloc=figure_dir)
    except:
        print(f"Error in along-shore visualisation {case=}")


### Figures - Contour
try:
    fv.vis_contour(data, t=t_moments, saveloc=figure_dir, variable="wl", xlims=[0, 1000], ylims=[0, 120])
    fv.vis_contour(data, t=t_moments, saveloc=figure_dir, variable="u", xlims=[0, 1000], ylims=[0, 120])
    fv.vis_contour(data, t=t_moments, saveloc=figure_dir, variable="v", xlims=[0, 1000], ylims=[0, 120])
except:
    print(f"Error in contour visualisation {case=}")


### Figures - Timeseries
try:
    fv.vis_timeseries(data, y=1e5*np.array([4, 8]), saveloc=figure_dir)
    fv.vis_timeseries(data, y=1e5*np.array([2, 4, 6, 8]), saveloc=figure_dir)
except:
    print(f"Error in timeseries visualisation {case=}")
    

### Spectra
## 1d
for _y in np.arange(1, 10):
    _y *= 1e5
    fv.vis_spectrum_1d(data, x=x_moment, y=_y, saveloc=figure_dir)
    fv.vis_spectrum_1d(data, x=x_moment+x_offset, y=_y, saveloc=figure_dir)

## 2d
fv.vis_spectrum_2d(data, x=x_moment, saveloc=figure_dir)
fv.vis_spectrum_2d(data, x=x_moment+x_offset, saveloc=figure_dir)


### Animation
if make_ani:
    anim.animation_contour(data, saveloc=figure_dir, xlims=[0, 1000])
    anim.animation_contour_uv(data, saveloc=figure_dir, xlims=[0, 1000])
    anim.animation_alongshore(data, saveloc=figure_dir, xlims=[0, 1000])
else:
    print("\nNot creating animations")


### End
with open(filename_log, "a+") as _file:
    _file.write(f"Created visualisations for case {case}\n")
        
if show_figs: 
    plt.show()
