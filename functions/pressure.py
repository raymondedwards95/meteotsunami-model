""" Functions to write time and space varying pressure fields for use with Delft3D-FM

Main functions:
    convert_to_xarray
    write_pressure
    plot_pressure
"""

import os
import sys
import time

import cmocean as cmo
import dask.array as da
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import numpy.typing as npt
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.utilities as fu
# fmt: on


def _find_step(
    x: npt.ArrayLike,
) -> Numeric:
    """ Find the step between two values in an array """
    dx = np.gradient(x)

    if not np.all(np.isclose(np.mean(dx), dx)):
        print(
            f"# Values are not equally spaced (mean: {np.mean(x)}; sd: {np.std(x)})")

    return np.mean(dx)


def convert_to_xarray(
    t: npt.ArrayLike,
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    p: npt.ArrayLike,
    savename: str = None,
    close: bool = False,
) -> xr.DataArray | None:
    """ Function to convert separate numpy arrays for t, x, y and p to a single DataArray for writing to files

    Input:
        `t`:        1d array with time since 1970-01-01 00:00:00 in seconds (shape = M)
        `x`:        1d array with x coordinates in km (shape = Nx)
        `y`:        1d array with y coordinates in km (shape = Ny)
        `p`:        3d array with pressure in Pa (shape = (M, Ny, Nx))

    Options:
        `savename`: saves data as an .nc file
        `close`:    closes data-array

    Output:
        `data`:     data-array containing data; if 'close=True' then data is None
    """
    print(f"# Converting pressure data to a data-array")
    t0 = time.perf_counter_ns()
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

    data = data.to_dataset(name="p")

    if savename is not None:
        if not savename.endswith(".nc"):
            savename += ".nc"
        data.to_netcdf(
            savename,
            encoding={"p": {
                "zlib": True,
                "complevel": 1,
                "least_significant_digit": 3
            }}
        )
        print(f"# Saved data-array as {savename}")

    t1 = time.perf_counter_ns()
    print(
        f"# Finished converting pressure data to a data-array in {(t1-t0)*1e-9:0.3f} seconds")

    if close:
        return
    return data["p"]


def filter_data(data: xr.DataArray) -> xr.DataArray:
    t0 = time.perf_counter_ns()

    data = data.round(2)
    ix = np.where(~ np.all(np.isclose(data, 0), axis=(0, 1)))[0]
    iy = np.where(~ np.all(np.isclose(data, 0), axis=(0, 2)))[0]
    data = data[:, :, ix][:, iy, :]

    t1 = time.perf_counter_ns()
    print(f"# Filtered data in {(t1-t0)*1e-9:0.3f} seconds")
    print(f"# New dimensions are {data.shape}")
    return data


def write_pressure(
    data: xr.DataArray,
    filename: str,
) -> None:
    """ Function to write a pressure field to a `pressure.amp` file for use with Delft3D-FM

    Input:
        `data`:     pressure field and coordinate data
        `filename`: name of .amp file
    """
    # prepare
    t0 = time.perf_counter_ns()
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if not filename.endswith(".amp"):
        filename += ".amp"
    print(f"# Writing pressure data to '{filename}'")

    # remove zero columns and rows
    data = filter_data(data)

    # prepare data
    x_num = data["x"].size
    y_num = data["y"].size
    t_num = data["t"].size

    x_min = data["x"].min().values
    y_min = data["y"].min().values

    dx = _find_step(data["x"].values)
    dy = _find_step(data["y"].values)

    x = data["x"].values
    y = data["y"].values
    t = data["t"].values
    p = data.chunk({"t": "auto", "x": -1, "y": "auto"}).data.compute()

    # header
    # fmt: off
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
    # fmt: on

    # write
    t_factor = np.min([t_num, 10])

    with open(filename, "w") as file:
        file.write(header)

        # loop over time
        for i in range(t_num):
            # fmt: off
            # progress
            if not (t_num-i-1) % (t_num // t_factor):
                print(f"# Step {i+1:4.0f} of {t_num:0.0f} ({(i+1)/t_num*100:0.1f}%)")

            file.write(
                f"TIME = {t[i]/3600.:0.06f} hours since 1970-01-01 00:00:00 +00:00\n".replace(".000000", ".0")
            )
            # fmt: on

            # put nothing if all values are close to reference
            # if ((i != 0) and (np.all(np.isclose(data.values[i, :, :].round(2), 0.)))):  # maybe set rtol and atol?
            # if (i > 0) and np.all(p[i,:,:].round(3) == p[i-1,:,:].round(3)):  # round to 3 numbers?
            # NOTE: it does seem it doesn't always work in d3dfm, so for now just write all data, even if it is all zeros
            if False:
                print(f"# \nSkip writing for hour {t[i]/3600.:0.3f} ({i=})")
                pass

            # write if there are values other than reference
            else:
                # loop over y (rows)
                for n in range(y_num):
                    # NOTE: negative y-index for p (i.e. -1*n) to fix coordinate-system
                    # replacement rules:
                    # 0.: values are rounded to two digits
                    # 1.: -0.00 -> 0.00 (remove minus)
                    # 2.: 0.00 -> 0 (remove .00)

                    p_ty = p[i, -1*n, :]
                    p_ty = p_ty.astype(str)
                    p_ty = np.char.replace(p_ty, "-0.0", "0.0")
                    p_ty = np.char.replace(p_ty, "0.0", "0")
                    s = " ".join(p_ty)

                    # write lines
                    file.write(s)
                    file.write("\n")

    # End
    t1 = time.perf_counter_ns()
    print(f"# Finished writing to '{filename}' in {(t1-t0)*1e-9:0.3f} seconds")

    return


def plot_pressure(
    data: xr.DataArray,
    filename: str = None,
    x_scales: Numeric = None,
    keep_open: bool = False,
) -> plt.Figure:
    """ Function to visualize pressure data

    Input:
        `data`:         pressure and coordinate data

    Options:
        `filename`:     name of figures
        `x_scales`:     lower and upper limit of x (should be a list of length 2)
        `keep_open`:    keep figures open after finishing
    """
    # prepare
    t0 = time.perf_counter_ns()
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if filename is None:
        filename = f"{os.path.dirname(os.path.realpath(__file__))}/tests/fig_pressure"
    if filename.endswith(".jpg"):
        filename.replace(".jpg", "")
    savename = f"{filename}_field"
    print(f"# Visualizing pressure field in '{savename}'")

    t_num = 5

    x = data["x"].values
    y = data["y"].values
    t = np.linspace(data["t"].min(), data["t"].max(), t_num)
    p = data.interp(t=t).compute().values

    p_max = np.ceil(np.max([p.max(), np.abs(p.min())]))
    p_min = -1. * p_max

    if x_scales is None:
        x_scales = [x.min() / 1000., x.max() / 1000.]

    # figure
    fig, ax = plt.subplots(1, t_num+1, sharey=True, squeeze=False,
                           gridspec_kw={"width_ratios": [4]*t_num + [1]})
    fig.set_size_inches(FIGSIZE_WIDE)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Pressure Disturbance", va="top", ha="left", x=0.01)
    ax = np.ravel(ax)

    im = [None] * t_num
    for i in range(t_num):
        idx = t.size * i // t_num
        im[i] = ax[i].contourf(
            x / 1000.,
            y / 1000.,
            p[i, :, :],
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
    cbar.set_ticks(np.linspace(np.floor(p.min()), np.ceil(p.max()), 11))

    fig.get_layout_engine().execute(fig)
    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI,
                pil_kwargs={"optimize": True, "compress_level": 9})
    if not keep_open:
        plt.close(fig)
    print(f"# Saved figure as '{savename}'")

    # End
    t1 = time.perf_counter_ns()
    print(f"# Finished visualising in {(t1-t0)*1e-9:0.3f} seconds")

    return fig


if __name__ == "__main__":
    # Define paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    pressure_dir = f"{script_dir}/tests/pres"
    pressure_file = f"{pressure_dir}/pressure"
    os.makedirs(pressure_dir, exist_ok=True)

    # Define function for computing 'pressure'
    def f(x, y, t):
        return (1. - da.exp(-t/(3.*3600.))) * da.sin(x * np.pi / 2000.) * da.sin(y * np.pi / 2000.)

    # Define grid
    x = np.linspace(-5000, 5000, 51)
    y = np.linspace(-10000, 10000, 51)
    t = np.linspace(0, 10, 15) * 3600.

    tt, yy, xx = da.meshgrid(t, y, x, indexing="ij")
    tt = tt.rechunk("auto")
    yy = yy.rechunk(tt.chunksize)
    xx = xx.rechunk(tt.chunksize)

    # Compute data
    p = f(xx, yy, tt)

    # Convert data
    convert_to_xarray(t, x, y, p, savename=pressure_file, close=True)
    del t, x, y, p

    # Read data from .nc-file
    data = xr.open_dataarray(f"{pressure_file}.nc", chunks="auto")

    # Write data to .amp-file
    write_pressure(data, filename=pressure_file)

    # Visualise data
    plot_pressure(data, filename=pressure_file)
