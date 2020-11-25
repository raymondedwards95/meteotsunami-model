""" Functions to process output files from Delft3D-FM into arrays """

import os

import numpy as np
import scipy.interpolate
import xarray as xr


def _convert_coordinate_to_regular_grid(coordinate, numsteps=None):
    if numsteps is None:
        numsteps = np.unique(coordinate).size
    return np.linspace(np.min(coordinate), np.max(coordinate), numsteps, retstep=True)


def _regrid_variable(var, x, y, x_grid, y_grid, index=None):
    assert var.ndim == 2

    if index is None:
        num_steps = var.shape[0]
    else:
        num_steps = 1

    var_regrid = np.zeros((num_steps, y_grid.size, x_grid.size), dtype=np.float)  # shape: time, y_grid, x_grid
    xy = np.vstack((x, y)).transpose()
    x_grid_mesh, y_grid_mesh = np.meshgrid(x_grid, y_grid)
    f_progress = np.max([num_steps // 10, 1])


    # loop over time
    for i in range(num_steps):
        # show progress
        if not (i+1) % f_progress:
            print(
                f"Step {i+1 : 3.0f}/{num_steps : 3.0f}")
        temp = scipy.interpolate.griddata(xy, var[i, :], (x_grid_mesh, y_grid_mesh), "linear")
        var_regrid[i, :, :] = temp

    return var_regrid


def extract_data(filename: str, savename: str = None):
    assert filename.endswith(".nc")
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
    data.t.attrs["long_name"] = "time"
    data.t.attrs["units"] = "-"

    ## Regrid variables
    print("Processing bathymetry")
    data["b"] = (("y", "x"), _regrid_variable(b, x, y, x_grid, y_grid, index=1)[0, :, :])
    data.b.attrs["long_name"] = "Water depth"
    data.b.attrs["units"] = "m"

    print("Processing water levels")
    data["wl"] = (("t", "y", "x"), _regrid_variable(wl, x, y, x_grid, y_grid))
    data.wl.attrs["long_name"] = "Water level"
    data.wl.attrs["units"] = "m"

    print("Processing zonal flow velocity")
    data["u"] = (("t", "y", "x"), _regrid_variable(u, x, y, x_grid, y_grid))
    data.u.attrs["long_name"] = "Zonal flow velocity"
    data.u.attrs["units"] = "m s-1"

    print("Processing meridional flow velocity")
    data["v"] = (("t", "y", "x"), _regrid_variable(v, x, y, x_grid, y_grid))
    data.v.attrs["long_name"] = "Meridional flow velocity"
    data.v.attrs["units"] = "m s-1"

    print("Processing atmospheric pressure")
    data["p"] = (("t", "y", "x"), _regrid_variable(p, x, y, x_grid, y_grid))
    data.p.attrs["long_name"] = "Atmospheric pressure near surface"
    data.p.attrs["units"] = "N m-2"

    print("Finished extracting data")

    if savename is not None:
        data.to_netcdf(savename)
        print(f"Data saved to '{savename}'")

    return data
