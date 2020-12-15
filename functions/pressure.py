""" Functions to write time and space varying pressure fields for use with Delft3D-FM """

import os

import numpy as np
import xarray as xr


def _find_step(x):
    """ Find the step between two values in an array """
    dx = np.gradient(x)

    if not np.all(np.isclose(np.mean(dx), dx)):
        print(f"WARNING: Values are not equally spaced (mean: {np.mean(x)}; sd: {np.std(x)})")

    return np.mean(dx)


def convert_to_xarray(t, x, y, p):
    """ Function to convert separate numpy arrays for t, x, y and p to a single DataArray for writing to files 
    
    Input:
        t:      1d array with time since 1970-01-01 00:00:00 in seconds (shape = M)
        x:      1d array with x coordinates in km (shape = Nx)
        y:      1d array with y coordinates in km (shape = Ny)
        p:      3d array with pressure in Pa (shape = (M, Ny, Nx))
    
    Output:
        data
    """
    data = xr.DataArray(p, dims=list("tyx"), coords={"x": x, "y": y, "t": t})

    data.attrs["long_name"] = "atmospheric pressure"
    data.attrs["units"] = "Pa"
    data.attrs["description"] = ""

    data.x.attrs["long_name"] = "x coordinate"
    data.x.attrs["units"] = "km"
    data.y.attrs["long_name"] = "y coordinate"
    data.y.attrs["units"] = "km"
    data.t.attrs["long_name"] = "seconds since 1970-01-01 00:00:00"
    data.t.attrs["units"] = "s"

    return data


def write_pressure(data, filename=None):
    """ Function to write a pressure field to a `pressure.amp` file for use with Delft3D-FM
    
    Input:
        data:   pressure field and coordinate data
    """
    ## prepare
    assert type(data) == xr.DataArray, "Input is not a DataArray"
    
    if filename is None:
        filename = os.path.dirname(os.path.realpath(__file__)) + "/pressure.amp"
    print(f"\nWriting pressure data to '{filename}'")

    x_num = data.x.values.size
    y_num = data.y.values.size
    t_num = data.t.values.size

    x_min = data.x.values.min()
    y_min = data.y.values.min()

    dx = _find_step(data.x.values)
    dy = _find_step(data.y.values)


    ## header
    header = f"""### START OF HEADER
FileVersion     = 1.03
filetype        = meteo_on_equidistant_grid
NODATA_value    = 0
n_cols          = {x_num:0.0f}
n_rows          = {y_num:0.0f}
grid_unit       = m
x_llcorner      = {x_min}
y_llcorner      = {y_min}
dx              = {dx}
dy              = {dy}
n_quantity      = 1
quantity1       = ap
unit1           = Pa
### END OF HEADER
"""


    ## write
    with open(filename, "w") as file:
        file.write(header)

        # loop over time
        for i in range(t_num):
            file.write(
                f"TIME = {data.t.values[i]/3600.:0.06f} hours since 1970-01-01 00:00:00 +00:00\n")

            # put nothing if all values are close to reference
            # if np.all(np.isclose(data.values[i, :, :].round(2), 0.)):  # maybe set rtol and atol?
            #     pass
            if False:
                pass

            # write if there are values other than reference
            else:

                # loop over y (rows)
                for n in range(y_num):

                    # loop over x (columns)
                    for m in range(x_num):

                        # write a value and neglect trailing zeros after decimal point
                        file.write(f"{data.values[i, n, m]:0.2f} ".replace(".00", ""))
                    file.write("\n")
    
    print(f"Finished writing to '{filename}'")


if __name__ == "__main__":
    print("")
