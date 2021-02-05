""" Functions to write bathymetry data in xyz- or xyb-format """

import os

import numpy as np
import xarray as xr


def _convert_to_xarray_1d(x, y, b):
    raise NotImplementedError()


def _convert_to_xarray_2d(x, y, b):
    return xr.DataArray(b, dims=list("yx"), coords={"x": x, "y": y})


def convert_to_xarray(x, y, b):
    """ Function to convert separate numpy arrays for x, y and b to a single DataArray for writing to files 
    
    Input:
        x:      1d array with x coordinates in km (shape = Nx)
        y:      1d array with y coordinates in km (shape = Ny)
        b:      1d or 2d array with bathymetry data in Pa (shape = (Ny * Nx) OR shape = (Ny, Nx))
    
    Output:
        data
    """
    if b.ndim == 1:
        data = _convert_to_xarray_1d(x, y, b)
    elif b.ndim == 2:
        data = _convert_to_xarray_2d(x, y, b)
    else:
        raise ValueError(f"{b.ndim} dimensions is not supported (max dim is 2")

    data.attrs["long_name"] = "bed level"
    data.attrs["units"] = "m"
    data.attrs["description"] = ""

    data.x.attrs["long_name"] = "x coordinate"
    data.x.attrs["units"] = "km"
    data.y.attrs["long_name"] = "y coordinate"
    data.y.attrs["units"] = "km"

    return data


def write_bathymetry(data, filename=None):
    """ Function to write bathymetry data to a `bethymetry.xyb` file for use with Delft3D-FM
    
    Input:
        data:   bathymetry and coordinate data
    """
    ## prepare
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if filename is None:
        filename = os.path.dirname(os.path.realpath(__file__)) + "/bathymetry"
    if not filename.endswith(".xyb"):
        filename += ".xyb"
    print(f"\nWriting bathymetry data to '{filename}'")

    x = data.x.values
    y = data.y.values
    b = data.values


    ## write
    with open(filename, "w") as file:
        file.write("")

        # loop over x
        for i in range(x.size):
            _x = x[i]

            # loop over y
            for j in range(y.size):

                # write a set of x y b on one line + remove trailing zeros
                file.write(f"{_x:0.2f} {y[j]:0.2f} {b[j, i]:0.2f}\n".replace(".00", ""))
    
    print(f"Finished writing to '{filename}'")


if __name__ == "__main__":
    import time

    @np.vectorize
    def function(x, y):
        return np.min([x / 500., 200.])

    t0 = time.perf_counter()
    x = np.linspace(0., 1e6, 5, dtype=np.float)
    y = np.linspace(-1e7, +1e7, 5, dtype=np.float)
    xx, yy = np.meshgrid(x, y)
    b = function(xx, yy)

    data = convert_to_xarray(x, y, b)
    write_bathymetry(data)

    t1 = time.perf_counter()
    print(f"Created bathymetry data in {t1 - t0 :0.2f} seconds.")
