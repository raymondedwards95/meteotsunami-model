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
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import numpy.typing as npt
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
# fmt: on


def _find_step(
    x: npt.ArrayLike,
) -> Numeric:
    """Find the step between two values in an array"""
    dx = np.gradient(x)

    if not np.all(np.isclose(np.mean(dx), dx)):
        print(f"# Values are not equally spaced (mean: {np.mean(x)}; sd: {np.std(x)})")

    return np.mean(dx)


def convert_to_xarray(
    t: npt.ArrayLike,
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    p: npt.ArrayLike,
    savename: str = None,
    close: bool = False,
) -> xr.DataArray | None:
    """Function to convert separate numpy arrays for t, x, y and p to a single DataArray for writing to files

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
            encoding={
                "p": {"zlib": True, "complevel": 1, "least_significant_digit": 3}
            },
        )
        print(f"# Saved data-array as {savename}")

    t1 = time.perf_counter_ns()
    print(
        f"# Finished converting pressure data to a data-array in {(t1-t0)*1e-9:0.3f} seconds"
    )

    if close:
        return
    return data["p"]


def filter_pressure(
    data: xr.DataArray,
    xmin: Numeric = None,
    xmax: Numeric = None,
    ymin: Numeric = None,
    ymax: Numeric = None,
    offset: Numeric = 1,
) -> xr.DataArray:
    """Rounds data and remove columns and rows that only contain zeros"""
    t0 = time.perf_counter_ns()
    print(f"# Start filtering data")
    shape = data.shape

    # Process options
    ixmin: int = -1
    ixmax: int = -1
    iymin: int = -1
    iymax: int = -1

    if xmin is not None:
        ixmin = np.searchsorted(data["x"].values, xmin)
    if xmax is not None:
        ixmax = np.searchsorted(data["x"].values, xmax)
    if ymin is not None:
        iymin = np.searchsorted(data["y"].values, ymin)
    if ymax is not None:
        iymax = np.searchsorted(data["y"].values, ymax)

    # Round
    data = data.round(2)

    # Find cols and rows where all values are close to 0
    ix = np.argwhere(~np.all(np.isclose(data, 0), axis=(0, 1)))
    iy = np.argwhere(~np.all(np.isclose(data, 0), axis=(0, 2)))

    # Apply options
    if xmin is not None:
        ixmin = np.min([ixmin, np.min(ix)])
    else:
        ixmin = np.min(ix)

    if xmax is not None:
        ixmax = np.max([ixmax, np.max(ix)])
    else:
        ixmax = np.max(ix)

    if ymin is not None:
        iymin = np.min([iymin, np.min(iy)])
    else:
        iymin = np.min(iy)

    if ymax is not None:
        iymax = np.max([iymax, np.max(iy)])
    else:
        iymax = np.max(iy)

    # Apply offset
    # Minimum index should be 0 or larger, certainly not negative
    # Maximum index does not have such restrictions, but limiting to size of original array is safe
    ixmin = np.max([0, ixmin - offset])
    ixmax = np.min([data["x"].size, ixmax + offset])
    iymin = np.max([0, iymin - offset])
    iymax = np.min([data["y"].size, iymax + offset])

    # Only remove cols and rows if they are on the outside
    slice_ix = slice(ixmin, ixmax + 1)
    slice_iy = slice(iymin, iymax + 1)

    # Apply slicing
    data = data[:, :, slice_ix][:, slice_iy, :]

    # End
    t1 = time.perf_counter_ns()
    print(f"# Old dimensions: {shape}")
    print(f"# New dimensions: {data.shape}")
    print(
        f"# ix = [{ixmin}, {ixmax}] -> x = [{data['x'].values.min():0.1f}, {data['x'].values.max():0.1f}]"
    )
    print(
        f"# iy = [{iymin}, {iymax}] -> x = [{data['y'].values.min():0.1f}, {data['y'].values.max():0.1f}]"
    )
    print(f"# Filtered data in {(t1-t0)*1e-9:0.3f} seconds")
    return data


def write_pressure(
    data: xr.DataArray,
    filename: str,
) -> None:
    """Function to write a pressure field to a `pressure.amp` file for use with Delft3D-FM

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
    data = filter_pressure(data)

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
    p = data.chunk({"t": "auto", "x": -1, "y": -1}).data.compute()

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
                f"TIME = {t[i]/3600.:0.6f} hours since 1970-01-01 00:00:00 +00:00\n".replace(".000000", ".0")
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

                    p_ty = p[i, -1 * n, :]
                    s = (
                        " ".join(f"{elem:0.02f}" for elem in p_ty)
                        .replace("-0.00", "0.00")
                        .replace("0.00", "0")
                    )

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
    x_min: Numeric = None,
    x_max: Numeric = None,
    y_min: Numeric = None,
    y_max: Numeric = None,
) -> tuple[plt.Figure, plt.Figure]:
    """Function to visualize pressure data

    Input:
        `data`:         pressure and coordinate data

    Options:
        `filename`:     name of figures
        `keep_open`:    keep figures open after finishing
        `x_min`:        lower limit of x
        `x_max`:        upper limit of x
        `y_min`:        lower limit of y
        `y_max`:        upper limit of y
    """
    # prepare
    t0 = time.perf_counter_ns()
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if filename is None:
        filename = f"{os.path.dirname(os.path.realpath(__file__))}/tests/fig_pressure"
    if filename.endswith(".jpg"):
        filename.replace(".jpg", "")

    t_num = 5

    # filter data
    if x_min is None:
        x_min = data["x"].values.min()
    data = filter_pressure(
        data,
        xmin=x_min,
        xmax=x_max,
        ymin=y_min,
        ymax=y_max,
    )

    # extract data
    x = data["x"].values
    y = data["y"].values
    t = np.linspace(data["t"].min(), data["t"].max(), t_num)
    p = data.interp(t=t).compute().values

    ix_single = data.argmax(["x", "y", "t"])["x"]

    p_max = np.ceil(np.max([p.max(), np.abs(p.min())]))

    if x_scales is None:
        x_scales = [x.min() / 1000.0, x.max() / 1000.0]

    cmap = cmo.cm.curl
    cbar_ticks = np.linspace(np.round(p.min()), np.round(p.max()), 11)

    # figure 1
    savename = f"{filename}_field"
    print(f"# Visualizing pressure field in '{savename}'")

    fig_1, ax_1 = plt.subplots(
        1,
        t_num + 1,
        sharey=True,
        squeeze=False,
        gridspec_kw={"width_ratios": [4] * t_num + [1]},
    )
    fig_1.set_size_inches(FIGSIZE_WIDE)
    fig_1.set_dpi(FIG_DPI)
    fig_1.set_layout_engine("compressed")
    fig_1.suptitle("Pressure Disturbance", va="top", ha="left", x=0.01)
    ax_1 = np.ravel(ax_1)

    im = [None] * t_num
    for i in range(t_num):
        idx = t.size * i // t_num
        im[i] = ax_1[i].pcolormesh(
            x / 1e3,
            y / 1e3,
            p[i, :, :],
            # levels=np.linspace(np.round(p.min()), np.round(p.max()), 100),
            cmap=cmap,
            vmin=-p_max,
            vmax=p_max,
            rasterized=True,
        )
        ax_1[i].set_title(f"$t = {t[idx]/3600.:0.0f}$h")
        ax_1[i].set_xlim(x_scales)  # make it automatic?

        # # rasterize contourf
        # for pathcoll in im[i].collections:
        #     pathcoll.set_rasterized(True)

    fig_1.supxlabel("$x$ [km]")
    fig_1.supylabel("$y$ [km]")

    # remove sharey from last subplot for colorbar
    ax_1[-1].get_shared_y_axes().remove(ax_1[-1])
    ax_1[-1].yaxis.major = matplotlib.axis.Ticker()
    ax_1[-1].yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    ax_1[-1].yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    cbar = fig_1.colorbar(im[-2], cax=ax_1[-1])
    cbar.set_label("Pressure Disturbance [Pa]")
    cbar.set_ticks(cbar_ticks)
    cbar.ax.set_ylim(cbar_ticks.min(), cbar_ticks.max())

    fig_1.get_layout_engine().execute(fig_1)
    save_figure(fig_1, savename)
    if not keep_open:
        plt.close(fig_1)

    # figure 2
    savename = f"{filename}_profile"
    print(f"# Visualizing pressure profile in '{savename}'")

    fig_2, ax_2 = plt.subplots(1, 1)
    fig_2.set_size_inches(FIGSIZE_WIDE)
    fig_2.set_dpi(FIG_DPI)
    fig_2.set_layout_engine("compressed")
    fig_2.suptitle("Pressure Disturbance", va="top", ha="left", x=0.01)

    for i in range(t_num):
        idx = t.size * i // t_num
        ax_2.plot(
            y / 1e3,
            p[i, :, ix_single],
            label=f"$t = {t[idx]/3600.:0.0f}$h",
        )
        ax_2.fill_between(
            y / 1e3,
            p[i, :, ix_single],
            alpha=0.1,
        )

    ax_2.axhline(
        color="black",
        linewidth=1,
        alpha=0.5,
    )

    ax_2.set_xlabel("$y$ [km]")
    ax_2.set_ylabel("Pressure Disturbance [Pa]")
    ax_2.grid()
    ax_2.legend()
    ax_2.set_xlim(y.min() / 1e3, y.max() / 1e3)
    ax_2.set_ylim(np.floor(p.min()), np.ceil(p.max()))

    fig_2.get_layout_engine().execute(fig_2)
    save_figure(fig_2, savename)
    if not keep_open:
        plt.close(fig_2)

    # End
    t1 = time.perf_counter_ns()
    print(f"# Finished visualising in {(t1-t0)*1e-9:0.3f} seconds")

    return fig_1, fig_2


if __name__ == "__main__":
    # Extra imports
    import dask.array as da

    # Define paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    pressure_dir = f"{PATH_TEST}/pressure"

    # Define function for computing 'pressure'
    def f(x, y, t):
        return (
            1.0
            * (1.0 - da.exp(-t / (3.0 * 3600.0)))
            * da.sin(x * np.pi / 2000.0)
            * da.sin(y * np.pi / 2000.0)
            * da.exp(-(y**2.0) / (3e3) ** 2.0)
        )

    # Define grid
    x = np.linspace(-5000, 5000, 101, dtype=np.float32)
    y = np.linspace(-15000, 15000, 200, dtype=np.float32)
    t = np.linspace(0, 10, 21, dtype=np.float32) * 3600.0

    tt, yy, xx = da.meshgrid(t, y, x, indexing="ij")
    tt = tt.rechunk("auto")
    yy = yy.rechunk(tt.chunksize)
    xx = xx.rechunk(tt.chunksize)

    # Compute data
    p = f(xx, yy, tt)

    # Convert data
    convert_to_xarray(t, x, y, p, savename=pressure_dir, close=True)
    del t, x, y, p

    # Read data from .nc-file
    data = xr.open_dataarray(f"{pressure_dir}.nc", chunks="auto")

    # Write data to .amp-file
    write_pressure(data, filename=pressure_dir)

    # Visualise data
    plot_pressure(data, filename=pressure_dir, y_min=-8e3)
