""" Functions to write time and space varying pressure fields for use with Delft3D-FM """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


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
        filename = os.path.dirname(os.path.realpath(__file__)) + "/tests/pressure"
    if not filename.endswith(".amp"):
        filename += ".amp"
    print(f"\nWriting pressure data to '{filename}'")

    x_num = data.x.values.size
    y_num = data.y.values.size
    t_num = data.t.values.size

    x_min = data.x.values.min()
    y_min = data.y.values.min()

    dx = _find_step(data.x.values)
    dy = _find_step(data.y.values)

    t = data.t.values
    p = data.values


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
    t_factor = np.min([t_num, 10])

    with open(filename, "w") as file:
        file.write(header)

        # loop over time
        for i in range(t_num):
            # progress
            if not (t_num-i-1) % (t_num // t_factor):
                print(f"Step {i+1:4.0f} of {t_num:0.0f} ({(i+1)/t_num*100:0.1f}%)")

            file.write(
                f"TIME = {t[i]/3600.:0.06f} hours since 1970-01-01 00:00:00 +00:00\n")

            # put nothing if all values are close to reference
            # if ((i != 0) and (np.all(np.isclose(data.values[i, :, :].round(2), 0.)))):  # maybe set rtol and atol?
            #     pass
            # if False:
            #     pass
            if (i > 0) and np.all(p[i,:,:].round(3) == p[i-1,:,:].round(3)):  # round to 3 numbers?
                print(f"\nSkip writing for hour {t[i]/3600.:0.3f} ({i=})")

            # write if there are values other than reference
            else:

                # loop over y (rows)
                for n in range(y_num):

                    # loop over x (columns)
                    for m in range(x_num):

                        # write a value and neglect trailing zeros after decimal point
                        # note: second index for p (i.e. -1*n) to fix coordinate-system
                        # replacement rules:
                        # 1.: -0 -> +0
                        # 2.: 0.00 -> 0
                        file.write(f"{p[i, -1*n, m]:0.2f} ".replace("-0.00", "0.00").replace(".00", ""))
                    file.write("\n")
    
    print(f"Finished writing to '{filename}'")


def plot_pressure(data, filename=None, x_scales=None):
    """ Function to visualize pressure data
    
    Input:
        data:       pressure and coordinate data
    
    Parameters:
        filename:   name of figures
        x_scales:   lower and upper limit of x (should be a list of length 2)
    """
    ## prepare
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if filename is None:
        filename = os.path.dirname(os.path.realpath(__file__)) + "/tests/fig_pressure"
    if filename.endswith(".jpg"):
        filename.replace(".jpg", "")
    savename = f"{filename}_field"
    print(f"\nVisualizing pressure field in '{savename}'")

    x = data.x.values
    y = data.y.values
    t = data.t.values
    p = data.values

    p_max = np.max([p.max(), np.abs(p.min())])
    p_min = -1. * p_max

    if x_scales is None:
        x_scales = [x.min() / 1000., x.max() / 1000.]

    ## figure
    t_num = 5
    fig, ax = plt.subplots(1, t_num+1, sharey=True, squeeze=False, gridspec_kw={"width_ratios": [4]*t_num + [1]})
    fig.set_size_inches(FIGSIZE_WIDE)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    ax = np.ravel(ax)

    im = [None] * t_num
    for i in range(t_num):
        idx = t.size * i // t_num
        im[i] = ax[i].contourf(
            x / 1000.,
            y / 1000.,
            p[idx, :, :],
            levels=31,
            vmin=p_min,
            vmax=p_max,
            cmap=cmo.cm.curl
        )
        ax[i].set_title(f"$t = {t[idx]/3600.:0.0f}$h")
        ax[i].set_xlim(x_scales)  # make it automatic?
    
    ax[t_num // 2].set_xlabel("$x$ [km]")
    ax[0].set_ylabel("$y$ [km]")

    ax[-1].get_shared_y_axes().remove(ax[-1])
    cbar = fig.colorbar(im[-2], cax=ax[-1])
    cbar.set_label("Pressure Disturbance [Pa]")

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI)
    print(f"Saved figure as '{savename}'")


if __name__ == "__main__":
    import time

    # function for computing 'pressure'
    def f(x, y, t):
        return np.min([1, t / 5. / 3600.]) * np.sin(x * np.pi) * np.sin(y * np.pi)

    t0 = time.perf_counter()

    # grid
    x = np.linspace(-5, 5, 51)
    y = np.linspace(-10, 10, 51)
    t = np.arange(0, 10) * 3600.
    p = np.zeros((t.size, y.size, x.size), dtype=np.float)

    # compute and convert data
    for n in range(t.size):
        for j in range(y.size):
            for i in range(x.size):
                p[n,j,i] = f(x[i], y[j], t[n])
    
    data = convert_to_xarray(t, x, y, p)
    del t, x, y, p

    # write data to .amp-file
    write_pressure(data)

    # evaluate performance
    t1 = time.perf_counter()
    print(f"Created test pressure data in {t1 - t0 :0.2f} seconds.")

    # visualise data
    plot_pressure(data)

    # evaluate performance
    t2 = time.perf_counter()
    print(f"Visualized test pressure data in {t2 - t1 :0.2f} seconds.")
