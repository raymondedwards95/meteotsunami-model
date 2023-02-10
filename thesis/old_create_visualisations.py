""" Script for making simple visualisations """
import argparse
import os
import sys
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.analysis as fa
import functions.animation as anim
import functions.regrid as fr
import functions.utilities as fu
import functions.visualisation as fv
# fmt: on


warnings.warn("This is an old script and should not be used!")


### Parameters
parser = argparse.ArgumentParser(
    description="Process model outputs and create visualisations for a specific case"
)
parser.add_argument("case", help="Case numbe", nargs="?", default="00", type=int)
parser.add_argument(
    "--reprocess",
    "-r",
    help="Force processing of original data.",
    default=False,
    type=bool,
)
parser.add_argument(
    "--show", "-s", help="Show figures after creation.", default=False, type=bool
)
parser.add_argument(
    "--animate", "-a", help="Create animation.", default=True, type=bool
)
parser.add_argument(
    "--delete-original-model-output",
    help="Delete original model output.",
    default=False,
    type=bool,
)
args = parser.parse_args()

case = f"{args.case:02}"
reprocess_data = bool(args.reprocess)
show_figs = bool(args.show)
make_ani = bool(args.animate)
delete_original_model_output = bool(args.delete_original_model_output)

if delete_original_model_output:
    warnings.warn("Original model output will be deleted!")
    print("Script will process data before deleting original model output!")
    reprocess_data = True


### Set paths
file_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = file_dir + "/output/"
output_dir = parent_dir + f"exp_{case}/"
figure_dir = file_dir + f"/visualisations_old/exp_{case}/"

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

# Delete original model output if asked
if delete_original_model_output:
    warnings.warn("Removing original model output!")
    time.sleep(2)
    os.remove(filename_output)
    print("Original model output is deleted!")
    with open(filename_log, "a+") as _file:
        _file.write(
            f"\tRemoved original model output for case {case}: {filename_output}\n"
        )

# Create folder for figures
os.makedirs(figure_dir, exist_ok=True)


### Get data
data = xr.open_dataset(filename_processed)

t_moments = np.array([8, 16, 24, 32, 40]) * 3600.0
x_moment = 1e4
x_offset = 5e4

y_list_0 = np.array([4, 8]) * 1e6
y_list_1 = np.array([2, 4, 6, 8]) * 1e6
y_moments = np.union1d(y_list_0, y_list_1)

waveperiods = np.zeros(y_moments.size)
wavelengths = np.zeros(t_moments.size)
wavespeeds = np.zeros((t_moments.size, y_moments.size))


### Compute stuff
with open(f"{figure_dir}/computed_parameters.txt", "w") as file:
    ## Waveperiod
    for i in range(y_moments.size):
        y_moment = y_moments[i]

        try:
            local_waveperiods = fa.compute_wave_periods(data, y_moment, x=x_moment)
            waveperiods[i] = np.nanmean(local_waveperiods)

            file.write(
                f"\n\nWave Period at y={y_moment/1000:0.0f} km and x={x_moment/1000:0.0f} km: {waveperiods[i]:0.1f} s\n"
            )
            for waveperiod in local_waveperiods:
                file.write(f"{waveperiod:0.1f}\t")

        except:
            print(f"*** Error in computing waveperiod {case=} and {y_moment=}")

    ## Wavelength
    for j in range(t_moments.size):
        t_moment = t_moments[j]

        try:
            local_wavelengths = fa.compute_wave_lengths(data, t_moment, x=x_moment)
            wavelengths[j] = np.nanmean(local_wavelengths)

            file.write(
                f"\n\nWave Lengths at t={t_moment/3600:0.1f} h and x={x_moment/1000:0.0f} km: {wavelengths[j]:0.1f} m\n"
            )
            for wavelength in local_wavelengths:
                file.write(f"{wavelength:0.1f}\t")

        except:
            print(f"*** Error in computing wavelength {case=} and {t_moment=}")

    ## Wavespeed
    file.write("\n\nWavespeeds:\n")
    for i in range(y_moments.size):
        for j in range(t_moments.size):
            try:
                wavespeeds[j, i] = wavelengths[j] / wavelengths[i]
                file.write(
                    f"y={y_moments[i]/1000:0.0f} km, t={t_moments[j]/3600:0.1f}h: {wavespeeds[j, i]:0.2f} m/s\n"
                )
            except:
                print(
                    f"*** Error in computing wavespeed for {case=}, {t_moments[j]=} and {y_moments[i]=}"
                )


### Figures
for j in range(t_moments.size):
    t_moment = t_moments[j]

    try:
        ## Find best parameters
        y_idx_max = fu.find_peaks_const_t(data, t_moment, x=x_moment, crests=True)
        y_idx_min = fu.find_peaks_const_t(data, t_moment, x=x_moment, crests=False)

        ## Cross-shore
        for i in range(np.min([y_idx_max.size, y_idx_min.size, 5])):
            fv.vis_crossshore(
                data, y=data["y"][y_idx_max[i]].values, t=t_moment, saveloc=figure_dir
            )
            fv.vis_crossshore(
                data, y=data["y"][y_idx_min[i]].values, t=t_moment, saveloc=figure_dir
            )
    except:
        print(f"*** Error in cross-shore visualisation {case=}")


### Figures - Contour
try:
    fv.vis_contour(
        data,
        t=t_moments,
        saveloc=figure_dir,
        variable="wl",
        xlims=[0, 1000],
        ylims=[0, 120],
    )
    fv.vis_contour(
        data,
        t=t_moments,
        saveloc=figure_dir,
        variable="u",
        xlims=[0, 1000],
        ylims=[0, 120],
    )
    fv.vis_contour(
        data,
        t=t_moments,
        saveloc=figure_dir,
        variable="v",
        xlims=[0, 1000],
        ylims=[0, 120],
    )
except:
    print(f"*** Error in contour visualisation {case=}")


### Figures - Alongshore
try:
    fv.vis_alongshore(data, x=x_moment, t=t_moments, saveloc=figure_dir)
    fv.vis_alongshore(data, x=x_moment + x_offset, t=t_moments, saveloc=figure_dir)
    for t_ in t_moments:
        fv.vis_alongshore(data, x=x_moment, t=t_, saveloc=figure_dir)
        fv.vis_alongshore(data, x=x_moment + x_offset, t=t_, saveloc=figure_dir)
except:
    print(f"*** Error in alongshore visualisation {case=}")


### Figures - Timeseries
try:
    fv.vis_timeseries(data, x=x_moment, y=y_list_0, saveloc=figure_dir)
    fv.vis_timeseries(data, x=x_moment, y=y_list_1, saveloc=figure_dir)
    fv.vis_timeseries(data, x=x_moment, y=y_moments, saveloc=figure_dir)
    fv.vis_timeseries(data, x=x_moment + x_offset, y=y_list_0, saveloc=figure_dir)
    fv.vis_timeseries(data, x=x_moment + x_offset, y=y_list_1, saveloc=figure_dir)
    fv.vis_timeseries(data, x=x_moment + x_offset, y=y_moments, saveloc=figure_dir)
except:
    print(f"*** Error in timeseries visualisation {case=}")


### Figures - Spectra 1d
try:
    for _y in y_moments:
        fv.vis_spectrum_1d(data, x=x_moment, y=_y, saveloc=figure_dir)
        fv.vis_spectrum_1d(data, x=x_moment + x_offset, y=_y, saveloc=figure_dir)
except:
    print(f"*** Error in spectrum-1d visualisation {case=}")

try:
    fv.vis_spectrum_1d(data, x=x_moment, y=y_list_0, saveloc=figure_dir)
    fv.vis_spectrum_1d(data, x=x_moment, y=y_list_1, saveloc=figure_dir)
    fv.vis_spectrum_1d(data, x=x_moment, y=y_moments, saveloc=figure_dir)
    fv.vis_spectrum_1d(data, x=x_moment + x_offset, y=y_list_0, saveloc=figure_dir)
    fv.vis_spectrum_1d(data, x=x_moment + x_offset, y=y_list_1, saveloc=figure_dir)
    fv.vis_spectrum_1d(data, x=x_moment + x_offset, y=y_moments, saveloc=figure_dir)
except:
    print(f"*** Error in spectrum-1d visualisation {case=}")

### Figures - Spectra 2d
try:
    fv.vis_spectrum_2d(data, x=x_moment, saveloc=figure_dir)
    fv.vis_spectrum_2d(data, x=+x_offset, saveloc=figure_dir)
except:
    print(f"*** Error in spectrum-2d visualisation {case=}")


### Figures - Cross-shore
try:
    fv.vis_crossshore(data, y=y_list_0, t=t_moments, saveloc=figure_dir)
    fv.vis_crossshore(data, y=y_list_1, t=t_moments, saveloc=figure_dir)
    fv.vis_crossshore(data, y=y_moments, t=t_moments, saveloc=figure_dir)

    for t_ in t_moments:
        fv.vis_crossshore(data, y=y_list_0, t=t_, saveloc=figure_dir)
        fv.vis_crossshore(data, y=y_list_1, t=t_, saveloc=figure_dir)
        fv.vis_crossshore(data, y=y_moments, t=t_, saveloc=figure_dir)

    for y_ in y_moments:
        fv.vis_crossshore(data, y=y_, t=t_moments, saveloc=figure_dir)

    for t_ in t_moments:
        for y_ in y_moments:
            fv.vis_crossshore(data, y=y_, t=t_, saveloc=figure_dir)
except:
    print(f"*** Error in cross-shore visualisation {case=}")


### Animation
if make_ani:
    anim.animation_alongshore(data, savedir=figure_dir, xlims=[0, 1000])
    anim.animation_contour(data, savedir=figure_dir, xlims=[0, 1000])
    anim.animation_contour_uv(data, savedir=figure_dir, xlims=[0, 1000])
else:
    print("\nNot creating animations")


### End
with open(filename_log, "a+") as _file:
    _file.write(f"Created visualisations for case {case}\n")

if show_figs:
    plt.show()
