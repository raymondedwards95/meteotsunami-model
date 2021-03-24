""" Functions to process output files from Delft3D-FM into arrays """

import os
import time
import warnings

import numpy as np
import scipy.interpolate
import xarray as xr


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
        grid_map = np.zeros((x.size, 3), dtype=np.int)  # slow method 1
    elif __regrid_method == 2:
        grid_map = np.zeros((x.size, 2), dtype=np.int)  # fast method 2
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

    if index is not None:
        warnings.warn("Warning: Parameter 'index' is not used here")
    
    # t
    t_size = var.shape[0]

    # x
    x_size = np.max(grid_mapping[:, 1]) + 1

    # y
    if __regrid_method == 1:
        y_size = np.max(grid_mapping[:, 2]) + 1  # slow method 1
    elif __regrid_method == 2:
        y_size = np.max(grid_mapping[:, 0]) + 1  # fast method 2
    else:
        raise NotImplementedError(__regrid_error)

    var_grid = np.full((t_size, y_size, x_size), np.nan, dtype=np.float)

    for k in range(t_size):
        # show progress
        # if not (k+1) % 1:
        #     print(f"Step {k+1 : 3.0f}/{t_size : 3.0f}")

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

    var_grid = np.zeros((num_steps, y_grid.size, x_grid.size), dtype=np.float)  # shape: time, y_grid, x_grid
    xy = np.vstack((x, y)).transpose()
    x_grid_mesh, y_grid_mesh = np.meshgrid(x_grid, y_grid)
    f_progress = np.max([num_steps // 10, 1])


    # loop over time
    for i in range(num_steps):
        # show progress
        if not (i+1) % f_progress:
            print(f"Step {i+1 : 3.0f}/{num_steps : 3.0f}")

        temp = scipy.interpolate.griddata(xy, var[i, :], (x_grid_mesh, y_grid_mesh), "linear")
        var_grid[i, :, :] = temp

    return var_grid


def extract_data(filename: str, savename: str = None):
    if not filename.endswith(".nc"):
        filename += ".nc"
    print("\nStart processing data")
    print(f"Reading file '{filename}'")

    ## Open data file and extract necessary coordinates and variables
    with xr.open_dataset(filename) as original_data:
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

    ## Regrid variables
    print("Processing bathymetry")
    t0 = time.perf_counter()
    # data["b"] = (("y", "x"), _regrid_variable_interpolate(b, x, y, x_grid, y_grid, index=1)[0, :, :])
    data["b"] = (("y", "x"), _regrid_variable_map(b, grid_map, index=1)[0, :, :])
    data.b.attrs["long_name"] = "Water depth"
    data.b.attrs["units"] = "m"
    t1 = time.perf_counter()
    print(f"  Used {t1 - t0:0.0f} seconds to process bathymetry data")

    print("Processing water levels")
    t0 = time.perf_counter()
    # data["wl"] = (("t", "y", "x"), _regrid_variable_interpolate(wl, x, y, x_grid, y_grid))
    data["wl"] = (("t", "y", "x"), _regrid_variable_map(wl, grid_map))
    data.wl.attrs["long_name"] = "Water level"
    data.wl.attrs["units"] = "m"
    t1 = time.perf_counter()
    print(f"  Used {t1 - t0:0.0f} seconds to process water level data")

    print("Processing zonal flow velocity")
    t0 = time.perf_counter()
    # data["u"] = (("t", "y", "x"), _regrid_variable_interpolate(u, x, y, x_grid, y_grid))
    data["u"] = (("t", "y", "x"), _regrid_variable_map(u, grid_map))
    data.u.attrs["long_name"] = "Zonal flow velocity"
    data.u.attrs["units"] = "m s-1"
    t1 = time.perf_counter()
    print(f"  Used {t1 - t0:0.0f} seconds to process zonal flow velocity data")

    print("Processing meridional flow velocity")
    t0 = time.perf_counter()
    # data["v"] = (("t", "y", "x"), _regrid_variable_interpolate(v, x, y, x_grid, y_grid))
    data["v"] = (("t", "y", "x"), _regrid_variable_map(v, grid_map))
    data.v.attrs["long_name"] = "Meridional flow velocity"
    data.v.attrs["units"] = "m s-1"
    t1 = time.perf_counter()
    print(f"  Used {t1 - t0:0.0f} seconds to process meridional flow velocity data")

    print("Processing atmospheric pressure")
    t0 = time.perf_counter()
    # data["p"] = (("t", "y", "x"), _regrid_variable_interpolate(p, x, y, x_grid, y_grid))
    data["p"] = (("t", "y", "x"), _regrid_variable_map(p, grid_map))
    data.p.attrs["long_name"] = "Atmospheric pressure near surface"
    data.p.attrs["units"] = "N m-2"
    t1 = time.perf_counter()
    print(f"  Used {t1 - t0:0.0f} seconds to process pressure data")

    print("\nFinished extracting data")

    if savename is not None:
        data.to_netcdf(savename)
        print(f"Data saved to '{savename}'")

    return data
