""" Functions to write time and space varying pressure fields for use with Delft3D-FM """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import xarray as xr

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


# make dir for tests
os.makedirs(f"{os.path.dirname(os.path.realpath(__file__))}/tests", exist_ok=True)


def _find_step(x):
    """ Find the step between two values in an array """
    dx = np.gradient(x)

    if not np.all(np.isclose(np.mean(dx), dx)):
        print(f"WARNING: Values are not equally spaced (mean: {np.mean(x)}; sd: {np.std(x)})")

    return np.mean(dx)


def convert_to_xarray(t, x, y, p, savename=None, close=False):
    """ Function to convert separate numpy arrays for t, x, y and p to a single DataArray for writing to files

    Input:
        t:          1d array with time since 1970-01-01 00:00:00 in seconds (shape = M)
        x:          1d array with x coordinates in km (shape = Nx)
        y:          1d array with y coordinates in km (shape = Ny)
        p:          3d array with pressure in Pa (shape = (M, Ny, Nx))

    Options:
        savename:   saves data as an .nc file

    Output:
        data:       xarray containing data; if 'close=True' then data is None
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

    if savename is not None:
        if not savename.endswith(".nc"):
            savename += ".nc"
        data.to_netcdf(
            savename,
            encoding={"__xarray_dataarray_variable__": {
                "zlib": True,
                "complevel": 1,
                "least_significant_digit": 3
            }}
        )

    if close:
        return
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
                f"TIME = {t[i]/3600.:0.06f} hours since 1970-01-01 00:00:00 +00:00\n".replace(".000000", ".0")
                )

            # put nothing if all values are close to reference
            # if ((i != 0) and (np.all(np.isclose(data.values[i, :, :].round(2), 0.)))):  # maybe set rtol and atol?
            # if (i > 0) and np.all(p[i,:,:].round(3) == p[i-1,:,:].round(3)):  # round to 3 numbers?
            # NOTE: it does seem it doesn't always work in d3dfm, so for now just write all data, even if it is all zeros
            if False:
                print(f"\nSkip writing for hour {t[i]/3600.:0.3f} ({i=})")
                pass

            # write if there are values other than reference
            else:

                # loop over y (rows)
                for n in range(y_num):

                    # loop over x (columns)
                    for m in range(x_num):
                        # write a value and neglect trailing zeros after decimal point
                        # NOTE: second index for p (i.e. -1*n) to fix coordinate-system
                        # replacement rules:
                        # 0.: values are rounded to two digits
                        # 1.: -0.00 -> 0.00 (remove minus)
                        # 2.: 0.00 -> 0 (remove .00)
                        file.write(f"{p[i, -1*n, m]:0.2f} ".replace("-0.00", "0.00").replace(".00", ""))
                    file.write("\n")

    print(f"Finished writing to '{filename}'")


def plot_pressure(data, filename=None, x_scales=None, keep_open=False):
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

    p_max = np.ceil(np.max([p.max(), np.abs(p.min())]))
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
            levels=41,
            vmin=p_min,
            vmax=p_max,
            cmap=cmo.cm.curl
        )
        ax[i].set_title(f"$t = {t[idx]/3600.:0.0f}$h")
        ax[i].set_xlim(x_scales)  # make it automatic?

    ax[t_num // 2].set_xlabel("$x$ [km]")
    ax[0].set_ylabel("$y$ [km]")

    # remove sharey from last subplot for colorbar
    ax[-1].get_shared_y_axes().remove(ax[-1])
    ax[-1].yaxis.major = matplotlib.axis.Ticker()
    ax[-1].yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    ax[-1].yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    cbar = fig.colorbar(im[-2], cax=ax[-1])
    cbar.set_label("Pressure Disturbance [Pa]")
    cbar.set_ticks(np.linspace(np.floor(p.min()), np.ceil(p.max()), 6))

    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    if not keep_open:
        plt.close(fig)
    print(f"Saved figure as '{savename}'")


if __name__ == "__main__":
    import time

    # function for computing 'pressure'
    def f(x, y, t):
        return np.min([1, t / 5. / 3600.]) * np.sin(x * np.pi / 2000.) * np.sin(y * np.pi / 2000.)

    t0 = time.perf_counter()

    # grid
    x = np.linspace(-5000, 5000, 51)
    y = np.linspace(-10000, 10000, 51)
    t = np.linspace(0, 10, 15) * 3600.
    p = np.zeros((t.size, y.size, x.size), dtype=float)

    # compute and convert data
    for n in range(t.size):
        for j in range(y.size):
            for i in range(x.size):
                p[n,j,i] = f(x[i], y[j], t[n])

    nc_file = f"{os.path.dirname(os.path.realpath(__file__))}/tests/pressure.nc"
    data = convert_to_xarray(t, x, y, p, savename=nc_file, close=True)
    del t, x, y, p, data

    # re-read data from .nc-file
    data = xr.open_dataarray(
        nc_file,
        chunks={
            "t": 1,
            "x": -1,
            "y": -1}
    )

    # write data to .amp-file
    write_pressure(data)
    data.close()

    # evaluate performance
    t1 = time.perf_counter()
    print(f"Created test pressure data in {t1 - t0 :0.2f} seconds.")

    # visualise data
    plot_pressure(data)

    # evaluate performance
    t2 = time.perf_counter()
    print(f"Visualized test pressure data in {t2 - t1 :0.2f} seconds.")

    #plt.show()
