""" Functions to process output files from Delft3D-FM into arrays """

import argparse
import os
import sys
import time
import warnings

import numpy as np
import scipy.interpolate
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


# set method for regridding
# 1 is the default; it should work
# 2 is faster, but can give wrong results
__regrid_method = 1
__regrid_error = f"Set '__regrid_method' in file 'regrid.py' to '1'! Now it is {__regrid_method}"


def _convert_coordinate_to_regular_grid(coordinate, numsteps=None):
    if numsteps is None:
        numsteps = np.unique(coordinate).size
    return np.linspace(np.min(coordinate), np.max(coordinate), numsteps, retstep=True)


def _create_grid_mapping(x, y, x_grid, y_grid):
    if type(x) is xr.DataArray:
        x = x.values
    if type(y) is xr.DataArray:
        y = y.values

    assert x.size == y.size
    # assert x.size == x_grid.size * y_grid.size

    if __regrid_method == 1:
        grid_map = np.zeros((x.size, 3), dtype=int)  # slow method 1
    elif __regrid_method == 2:
        grid_map = np.zeros((x.size, 2), dtype=int)  # fast method 2
    else:
        raise NotImplementedError(__regrid_error)

    for n in range(x.size):
        i = np.argmin(np.abs(x[n] - x_grid))
        j = np.argmin(np.abs(y[n] - y_grid))

        if __regrid_method == 1:
            grid_map[n, :] = [n, i, j]  # slow method 1
        elif __regrid_method == 2:
            grid_map[n, :] = [j, i]  # fast method 2
        else:
            raise NotImplementedError(__regrid_error)

    return grid_map


def _regrid_variable_map(var, grid_mapping, index=None):
    if type(var) is xr.DataArray:
        var = var.values
    assert var.ndim == 2
    assert grid_mapping.ndim == 2

    # t
    t_size = var.shape[0]
    if index is not None:
        warnings.warn("Warning: Parameter 'index' is not used properly in this function. Check source for details!")
        t_size = 1

    # x
    x_size = np.max(grid_mapping[:, 1]) + 1

    # y
    if __regrid_method == 1:
        y_size = np.max(grid_mapping[:, 2]) + 1  # slow method 1
    elif __regrid_method == 2:
        y_size = np.max(grid_mapping[:, 0]) + 1  # fast method 2
    else:
        raise NotImplementedError(__regrid_error)

    var_grid = np.full((t_size, y_size, x_size), np.nan, dtype=float)

    progress_factor = np.min([5, t_size])

    for k in range(t_size):
        # show progress
        if not (t_size-k-1) % (t_size // progress_factor):
            print(f"Step {k:4.0f} of {t_size:0.0f} ({(k+1)/t_size*100:0.1f}%)")

        if __regrid_method == 1:
            for m in range(grid_mapping.shape[0]):  # slow method 1
                n, i, j = grid_mapping[m, :]
                var_grid[k, j, i] = var[k, n]

        elif __regrid_method == 2:
            var_grid[k] = var[k, :].reshape(y_size, x_size)[tuple(grid_mapping.T)].reshape(y_size, x_size)  # fast method 2

        else:
            raise NotImplementedError(__regrid_error)

    return var_grid


def _regrid_variable_interpolate(var, x, y, x_grid, y_grid, index=None):
    warnings.warn("Warning: Using interpolation as regridding tool!")

    if type(x) is xr.DataArray:
        x = x.values
    if type(y) is xr.DataArray:
        y = y.values
    if type(var) is xr.DataArray:
        var = var.values
    assert var.ndim == 2

    if index is None:
        num_steps = var.shape[0]
    else:
        num_steps = 1

    var_grid = np.zeros((num_steps, y_grid.size, x_grid.size), dtype=float)  # shape: time, y_grid, x_grid
    xy = np.vstack((x, y)).transpose()
    x_grid_mesh, y_grid_mesh = np.meshgrid(x_grid, y_grid)

    progress_factor = np.min([5, num_steps])


    # loop over time
    for i in range(num_steps):
        # show progress
        if not (num_steps-i-1) % (num_steps // progress_factor):
            print(f"Step {i:4.0f} of {num_steps:0.0f} ({(i+1)/num_steps*100:0.1f}%)")

        temp = scipy.interpolate.griddata(xy, var[i, :], (x_grid_mesh, y_grid_mesh), "linear")
        var_grid[i, :, :] = temp

    return var_grid


def extract_data(filename: str, savename: str = None):
    t_start = time.perf_counter()

    if not filename.endswith(".nc"):
        filename += ".nc"
    print("\nStart processing data")
    print(f"Reading file '{filename}'")

    ## Open data file and extract necessary coordinates and variables
    with xr.open_dataset(filename) as original_data:
        x = original_data["mesh2d_face_x"]
        y = original_data["mesh2d_face_y"]
        t = original_data["time"]

        # b = original_data["mesh2d_waterdepth"]  # moved to part with '## Regrid variables'
        # wl = original_data["mesh2d_s1"]
        # u = original_data["mesh2d_ucx"]
        # v = original_data["mesh2d_ucy"]
        # p = original_data["mesh2d_Patm"]

    ## Make regular grid
    x_grid = _convert_coordinate_to_regular_grid(x)[0]
    y_grid = _convert_coordinate_to_regular_grid(y)[0]

    data = xr.Dataset(coords={"t": t.values, "x": x_grid, "y": y_grid})
    data.x.attrs["long_name"] = "x coordinate"
    data.x.attrs["units"] = "m"
    data.y.attrs["long_name"] = "y coordinate"
    data.y.attrs["units"] = "m"
    # data.t.attrs["long_name"] = "time"
    # data.t.attrs["units"] = "s"

    grid_map = _create_grid_mapping(x, y, x_grid, y_grid)

    ## Regrid variables
    print("\nProcessing bathymetry")
    t0 = time.perf_counter()
    with xr.open_dataset(filename) as original_data:
        b = original_data["mesh2d_waterdepth"]
    data["b"] = (("y", "x"), _regrid_variable_map(b, grid_map, index=1)[0, :, :])
    data.b.attrs["long_name"] = "Water depth"
    data.b.attrs["units"] = "m"
    del b
    t1 = time.perf_counter()
    print(f"Used {t1 - t0:0.0f} seconds to process bathymetry data")

    print("\nProcessing water levels")
    t0 = time.perf_counter()
    with xr.open_dataset(filename) as original_data:
        wl = original_data["mesh2d_s1"]
    data["wl"] = (("t", "y", "x"), _regrid_variable_map(wl, grid_map))
    data.wl.attrs["long_name"] = "Water level"
    data.wl.attrs["units"] = "m"
    del wl
    t1 = time.perf_counter()
    print(f"Used {t1 - t0:0.0f} seconds to process water level data")

    print("\nProcessing zonal flow velocity")
    t0 = time.perf_counter()
    with xr.open_dataset(filename) as original_data:
        u = original_data["mesh2d_ucx"]
    data["u"] = (("t", "y", "x"), _regrid_variable_map(u, grid_map))
    data.u.attrs["long_name"] = "Zonal flow velocity"
    data.u.attrs["units"] = "m s-1"
    del u
    t1 = time.perf_counter()
    print(f"Used {t1 - t0:0.0f} seconds to process zonal flow velocity data")

    print("\nProcessing meridional flow velocity")
    t0 = time.perf_counter()
    with xr.open_dataset(filename) as original_data:
        v = original_data["mesh2d_ucy"]
    data["v"] = (("t", "y", "x"), _regrid_variable_map(v, grid_map))
    data.v.attrs["long_name"] = "Meridional flow velocity"
    data.v.attrs["units"] = "m s-1"
    del v
    t1 = time.perf_counter()
    print(f"Used {t1 - t0:0.0f} seconds to process meridional flow velocity data")

    print("\nProcessing atmospheric pressure")
    t0 = time.perf_counter()
    with xr.open_dataset(filename) as original_data:
        p = original_data["mesh2d_Patm"]
    data["p"] = (("t", "y", "x"), _regrid_variable_map(p, grid_map))
    data.p.attrs["long_name"] = "Atmospheric pressure near surface"
    data.p.attrs["units"] = "N m-2"
    del p
    t1 = time.perf_counter()
    print(f"Used {t1 - t0:0.0f} seconds to process pressure data")

    print("\nFinished extracting data")

    if savename is not None:
        t0 = time.perf_counter()
        encoding = {var: {"zlib": True, "complevel": 1, "least_significant_digit": 9} for var in data}
        data.to_netcdf(savename, encoding=encoding)
        t1 = time.perf_counter()
        print(f"Data saved to '{savename}' in {t1 - t0:0.0f} seconds")

    t_end = time.perf_counter()
    t_total = t_end - t_start
    print(f"Finished regridding in {t_total:0.0f} seconds ({t_total/60:0.0f} minutes)")

    return data


if __name__ == "__main__":
    ### Options
    parser = argparse.ArgumentParser(
        description="Process and regrid model output"
    )
    parser.add_argument(
        "input",
        help="Name of the model output file (default file name is 'FlowFM_map.nc')",
        type=str
    )
    parser.add_argument(
        "output",
        help="Name of the (new) file to write the processed data",
        type=str
    )
    parser.add_argument(
        "--delete-original-model-output",
        help="Delete original model output",
        default=False,
        type=bool
    )
    args = parser.parse_args()

    ### Filenames
    filename_original = str(args.input)
    if not filename_original.endswith(".nc"):
        filename_original += ".nc"
    print(f"Input is '{filename_original}'")

    filename_processed = str(args.output)
    if not filename_processed.endswith(".nc"):
        filename_processed += ".nc"
    print(f"Output is '{filename_processed}'")

    delete_original_model_output = bool(args.delete_original_model_output)
    if delete_original_model_output:
        warnings.warn(f"Original model output will be deleted! {filename_original}")
        time.sleep(2)

    ### Convert data
    print("Processing data")
    time.sleep(2)
    extract_data(
        filename=filename_original,
        savename=filename_processed
    )
    print("Finished processing data")

    ### Clean up data
    if delete_original_model_output:
        warnings.warn(f"Original model output will be deleted! {filename_original}")
        time.sleep(2)
        os.remove(filename_original)
        print(f"Original model output is deleted! {filename_original}")
