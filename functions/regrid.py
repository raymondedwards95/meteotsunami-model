""" Functions to process output files from Delft3D-FM into structured arrays """

import argparse
import os
import sys
import time
import warnings
from typing import Union, Tuple

import numpy as np
import numpy.typing as npt
import scipy.interpolate
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


# set method for regridding
# 1 is the default; it should work
# 2 is faster, but sometimes it gives wrong results
__regrid_method = 1
__regrid_error = f"Set '__regrid_method' in file 'regrid.py' to '1'! Now it is {__regrid_method}"


def _convert_coordinate_to_regular_grid(coordinate: npt.ArrayLike, numsteps: Integer=None) -> Tuple[np.ndarray, Numeric]:
    """ Creates an equally spaced array of coordinates

    Input:
        coordinate: array with coordinates

    Options:
        numsteps:   number of coordinates

    Output:
        grid:       array of gridded coordinates
        step:       value of spacing in coordinates
    """
    if numsteps is None:
        numsteps = np.unique(coordinate).size
    grid, step = np.linspace(np.min(coordinate), np.max(coordinate), numsteps, retstep=True)
    return (grid, step)


def _create_grid_mapping(x: Union[xr.DataArray, npt.ArrayLike], y: Union[xr.DataArray, npt.ArrayLike], x_grid: npt.ArrayLike, y_grid: npt.ArrayLike) -> np.ndarray:
    """ Find a mapping to convert data in unstructured data to gridded data

    Input:
        x:          1-d array of unstructured x-coordinates
        y:          1-d array of corresponding y-coordinates
        x_grid:     1-d array of x-coordinates of a grid
        y_grid:     1-d array of y-coordinates of a grid

    Output:
        grid_map:   map that maps the unstructered coordinates to gridded coordinates
    """
    # Check inputs
    if type(x) is xr.DataArray:
        x = x.values
    if type(y) is xr.DataArray:
        y = y.values

    assert np.size(x) == np.size(y)
    # assert np.size(x) == np.size(x_grid) * np.size(y_grid)  # check would be nice, but not necessary

    # Prepare arrays
    if __regrid_method == 1:
        grid_map = np.zeros((x.size, 3), dtype=int)  # slow method 1
    elif __regrid_method == 2:
        grid_map = np.zeros((x.size, 2), dtype=int)  # fast method 2
    else:
        raise NotImplementedError(__regrid_error)

    # Find mapping, loop over all x-y pairs
    for n in range(x.size):
        i = np.argmin(np.abs(x[n] - x_grid))
        j = np.argmin(np.abs(y[n] - y_grid))

        if __regrid_method == 1:
            grid_map[n, :] = [n, i, j]  # slow method 1
        elif __regrid_method == 2:
            grid_map[n, :] = [j, i]  # fast method 2
        else:
            raise NotImplementedError(__regrid_error)

    # End
    return grid_map


def _regrid_variable_map(var: Union[xr.DataArray, npt.ArrayLike], grid_map: npt.ArrayLike, index: Integer=None) -> np.ndarray:
    """ Regrids unstructured data using a pre-defined mapping

    Input:
        var:        unstructured data
        grid_map:   map that maps the unstructered coordinates to gridded coordinates

    Options:
        index:      only regrid a specific slice

    Output:
        var_grid:   regridded data
    """
    # Check inputs
    if type(var) is xr.DataArray:
        var = var.values
    assert np.ndim(var) == 2
    assert np.ndim(grid_map) == 2


    # Prepare t
    t_size = var.shape[0]
    if index is not None:
        warnings.warn("Warning: Parameter 'index' is used!")
        t_size = 1
    progress_factor = np.min([5, t_size])

    # Prepare x
    x_size = np.max(grid_map[:, 1]) + 1

    # Prepare y
    if __regrid_method == 1:
        y_size = np.max(grid_map[:, 2]) + 1  # slow method 1
    elif __regrid_method == 2:
        y_size = np.max(grid_map[:, 0]) + 1  # fast method 2
    else:
        raise NotImplementedError(__regrid_error)

    # Prepare grid
    var_grid = np.full((t_size, y_size, x_size), np.nan, dtype=float)

    # Regrid loop over time
    for k in range(t_size):
        # Show progress
        if not (t_size-k-1) % (t_size // progress_factor):
            print(f"# Step {k:4.0f} of {t_size:0.0f} ({(k+1)/t_size*100:0.1f}%)")

        # Regrid
        if __regrid_method == 1:
            for m in range(grid_map.shape[0]):  # slow method 1
                n, i, j = grid_map[m, :]
                var_grid[k, j, i] = var[k, n]

        elif __regrid_method == 2:
            var_grid[k] = var[k, :].reshape(y_size, x_size)[tuple(grid_map.T)].reshape(y_size, x_size)  # fast method 2

        else:
            raise NotImplementedError(__regrid_error)

    # End
    return var_grid


def _regrid_variable_interpolate(var: Union[xr.DataArray, npt.ArrayLike], x: Union[xr.DataArray, npt.ArrayLike], y: Union[xr.DataArray, npt.ArrayLike], x_grid: npt.ArrayLike, y_grid: npt.ArrayLike, index: Integer=None) -> np.ndarray:
    """ Regrids unstructured data using interpolation

    Input:
        var:        unstructured data
        x:          1-d array of unstructured x-coordinates
        y:          1-d array of corresponding y-coordinates
        x_grid:     1-d array of x-coordinates of a grid
        y_grid:     1-d array of y-coordinates of a grid

    Options:
        index:      only regrid a specific slice

    Output:
        var_grid:   regridded data
    """
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
            print(f"# Step {i:4.0f} of {num_steps:0.0f} ({(i+1)/num_steps*100:0.1f}%)")

        temp = scipy.interpolate.griddata(xy, var[i, :], (x_grid_mesh, y_grid_mesh), "linear")
        var_grid[i, :, :] = temp

    return var_grid


def extract_data(filename: str, savename: str, close: bool=False) -> xr.DataSet | None:
    """ Extracts unstructured data from the output of D3D-FLOW-FM and convert it to a structured dataset

    Input:
        filename    name of the output file from D3D-FLOW-FM
        savename    name of the new dataset

    Options:
        close:      close dataset after finishing

    Output:
        data:       data array with structured data; if `close=False` then data is `None`
    """
    t0 = time.perf_counter_ns()

    if not filename.endswith(".nc"):
        filename += ".nc"
    if not savename.endswith(".nc"):
        savename += ".nc"

    print(f"#")
    print(f"# Start processing data")
    print(f"# Reading file '{filename}'")

    ## Open data file and extract necessary coordinates and variables
    with xr.open_dataset(filename, chunks="auto") as original_data:
        x = original_data["mesh2d_face_x"]
        y = original_data["mesh2d_face_y"]
        t = original_data["time"]

        b = original_data["mesh2d_waterdepth"]
        wl = original_data["mesh2d_s1"]
        u = original_data["mesh2d_ucx"]
        v = original_data["mesh2d_ucy"]
        p = original_data["mesh2d_Patm"]

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

    ## Regrid bathymetry
    print(f"#")
    print(f"# Processing bathymetry")
    ta = time.perf_counter_ns()
    data["b"] = (("y", "x"), _regrid_variable_map(b, grid_map, index=1)[0, :, :])
    data.b.attrs["long_name"] = "Water depth"
    data.b.attrs["units"] = "m"
    del b
    tb = time.perf_counter_ns()
    print(f"# Used {(tb-ta)*1e-9:0.3f} seconds to process bathymetry data")

    ## Regrid water levels
    print(f"#")
    print(f"# Processing water levels")
    ta = time.perf_counter_ns()
    data["wl"] = (("t", "y", "x"), _regrid_variable_map(wl, grid_map))
    data.wl.attrs["long_name"] = "Water level"
    data.wl.attrs["units"] = "m"
    del wl
    tb = time.perf_counter_ns()
    print(f"# Used {(tb-ta)*1e-9:0.3f} seconds to process water level data")

    ## Regrid zonal flow velocity
    print(f"#")
    print(f"# Processing zonal flow velocity")
    ta = time.perf_counter_ns()
    data["u"] = (("t", "y", "x"), _regrid_variable_map(u, grid_map))
    data.u.attrs["long_name"] = "Zonal flow velocity"
    data.u.attrs["units"] = "m s-1"
    del u
    tb = time.perf_counter_ns()
    print(f"# Used {(tb-ta)*1e-9:0.3f} seconds to process zonal flow velocity data")

    ## Regrid meridional flow velocity
    print(f"#")
    print(f"# Processing meridional flow velocity")
    ta = time.perf_counter_ns()
    data["v"] = (("t", "y", "x"), _regrid_variable_map(v, grid_map))
    data.v.attrs["long_name"] = "Meridional flow velocity"
    data.v.attrs["units"] = "m s-1"
    del v
    tb = time.perf_counter_ns()
    print(f"# Used {(tb-ta)*1e-9:0.3f} seconds to process meridional flow velocity data")

    ## Regrid atmospheric pressure
    print(f"#")
    print(f"# Processing atmospheric pressure")
    ta = time.perf_counter_ns()
    data["p"] = (("t", "y", "x"), _regrid_variable_map(p, grid_map))
    data.p.attrs["long_name"] = "Atmospheric pressure near surface"
    data.p.attrs["units"] = "N m-2"
    del p
    tb = time.perf_counter_ns()
    print(f"# Used {(tb-ta)*1e-9:0.3f} seconds to process pressure data")

    ## End of regrid
    print(f"#")
    print(f"# Finished extracting data")

    ## Save data
    ta = time.perf_counter_ns()
    encoding = {var: {"zlib": True, "complevel": 1, "least_significant_digit": 6} for var in data}
    data.to_netcdf(savename, encoding=encoding)
    data.close()
    tb = time.perf_counter_ns()
    print(f"# Data saved to '{savename}' in {(t1-t0)*1e-9:0.3f} seconds")

    ## End
    t1 = time.perf_counter_ns()
    t_total = (t1 - t0) * 1e-9
    print(f"#")
    print(f"# Finished regridding in {t_total:0.3f} seconds ({t_total/60:0.1f} minutes)")
    print(f"#")

    if close:
        return
    return xr.open_dataset(savename)


if __name__ == "__main__":
    ### Start
    print(f"\nStart time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

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
    print(f"Processing data")
    time.sleep(2)
    extract_data(
        filename=filename_original,
        savename=filename_processed
    )
    print(f"Finished processing data")

    ### Clean up data
    if delete_original_model_output:
        warnings.warn(f"Original model output will be deleted! {filename_original}")
        time.sleep(2)
        os.remove(filename_original)
        print(f"Original model output is deleted! {filename_original}")

    ### End
    print(f"\nFinish time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
