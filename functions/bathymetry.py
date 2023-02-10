""" Functions to write bathymetry data in xyz- or xyb-format

Main functions:
    convert_to_xarray
    write_bathymetry
    plot_bathymetry
"""

import os
import sys
import time

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
# fmt: on


def _convert_to_xarray_1d(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    b: npt.ArrayLike,
) -> xr.DataArray:
    if np.ndim(x) != 1:
        raise ValueError(f"'x' should be 1-dimensional, instead of {np.ndim(x)}")
    if np.ndim(y) != 1:
        raise ValueError(f"'y' should be 1-dimensional, instead of {np.ndim(y)}")
    if np.ndim(b) != 1:
        raise ValueError(f"'b' should be 1-dimensional, instead of {np.ndim(b)}")
    raise NotImplementedError()


def _convert_to_xarray_2d(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    b: npt.ArrayLike,
) -> xr.DataArray:
    if np.ndim(x) != 1:
        raise ValueError(f"'x' should be 1-dimensional, instead of {np.ndim(x)}")
    if np.ndim(y) != 1:
        raise ValueError(f"'y' should be 1-dimensional, instead of {np.ndim(y)}")
    if np.ndim(b) != 2:
        raise ValueError(f"'b' should be 2-dimensional, instead of {np.ndim(b)}")
    return xr.DataArray(b, dims=list("yx"), coords={"x": x, "y": y})


def convert_to_xarray(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    b: npt.ArrayLike,
    savename: str = None,
    close: bool = False,
) -> xr.DataArray | None:
    """Function to convert separate numpy arrays for x, y and b to a single DataArray for writing to files

    Input:
        `x`:        1d array with x coordinates in km (shape = Nx)
        `y`:        1d array with y coordinates in km (shape = Ny)
        `b`:        1d or 2d array with bathymetry data in m (shape = (Ny * Nx) OR shape = (Ny, Nx))

    Options:
        `savename`: saves data as an .nc file
        `close`:    closes data, if 'close=False', then data will be returned

    Output:
        `data`:     data-array containing data; if 'close=True' then data is None
    """
    print(f"# Converting bathymetry data to a data-array")
    t0 = time.perf_counter_ns()

    if b.ndim == 1:
        data = _convert_to_xarray_1d(x, y, b)
    elif b.ndim == 2:
        data = _convert_to_xarray_2d(x, y, b)
    else:
        raise ValueError(f"{b.ndim} dimensions is not supported for 'b' (max dim is 2)")

    data.attrs["long_name"] = "bed level"
    data.attrs["units"] = "m"
    data.attrs["description"] = ""

    data.x.attrs["long_name"] = "x coordinate"
    data.x.attrs["units"] = "km"
    data.y.attrs["long_name"] = "y coordinate"
    data.y.attrs["units"] = "km"

    if savename is not None:
        if not savename.endswith(".nc"):
            savename += ".nc"
        data.to_netcdf(
            savename,
            encoding={
                "__xarray_dataarray_variable__": {
                    "zlib": True,
                    "complevel": 1,
                    "least_significant_digit": 3,
                }
            },
        )
        print(f"# Saved data-array as {savename}")

    t1 = time.perf_counter_ns()
    print(
        f"# Finished converting bathymetry data to a data-array in {(t1-t0)*1e-9:0.3f} seconds"
    )

    if close:
        return
    return data


def write_bathymetry(
    data: xr.DataArray,
    filename: str,
) -> None:
    """Function to write bathymetry data to a `bathymetry.xyb` file for use with Delft3D-FM

    Input:
        `data`:     bathymetry and coordinate data
        `filename`: .xyb file name
    """
    # prepare
    t0 = time.perf_counter_ns()
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if not filename.endswith(".xyb"):
        filename += ".xyb"
    print(f"# Writing bathymetry data to '{filename}'")

    x = data.x.values
    y = data.y.values
    b = data.values

    # convert data to plaintext
    textdata = str()

    for i in range(x.size):
        # loop over x
        _x = x[i]

        for j in range(y.size):
            # loop over y

            # write a set of x y b on one line + remove trailing zeros
            textdata += f"{_x:0.2f} {y[j]:0.2f} {b[j, i]:0.2f}\n".replace(".00", "")

    # write
    with open(filename, "w") as file:
        file.write(textdata)

    t1 = time.perf_counter_ns()
    print(f"# Finished writing to '{filename}' in {(t1-t0)*1e-9:0.3f} seconds")


def plot_bathymetry(
    data: xr.DataArray,
    filename: str = None,
    xmax: Numeric = None,
    keep_open: bool = False,
    half_width: bool = False,
    scale: str = "Mm",
) -> tuple[plt.Figure, plt.Figure]:
    """Function to visualize bathymetry data

    Input:
        `data`:         bathymetry and coordinate data

    Options:
        `filename`:     name of figures
        `xmax`:         upper limit of x
        `keep_open`:    keep figures open after finishing
        `half_width`:   make figures for use in columns, i.e. side-by-side
        `scale`:        scale of plots ('m', 'km' or 'Mm')
    """
    # Prepare
    t0 = time.perf_counter_ns()
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if filename is None:
        script_path = os.path.dirname(os.path.realpath(__file__))
        filename = f"{script_path}/tests/fig_bathymetry"
    if filename.endswith(".jpg"):
        filename.replace(".jpg", "")

    figsize = FIGSIZE_NORMAL
    if half_width:
        figsize = FIGSIZE_SMALL
        filename += "_small"

    unit: str
    scale_factor: float

    match scale:
        case "m":
            unit = "\si{\meter}"
            scale_factor = 1e0
        case "km":
            unit = "\si{\kilo\meter}"
            scale_factor = 1e3
        case "Mm":
            unit = "\si{\mega\meter}"
            scale_factor = 1e6
        case _:
            raise ValueError(
                f"Scale should be either 'm', 'km', or 'Mm', instead of '{scale}'"
            )

    print(f"# Visualizing bathymetry in '{filename}'")

    # Extract data
    x = data.x.values
    y = data.y.values
    b = data.values

    b_max = np.ceil(b.max())
    b_min = np.floor(b.min())

    # fix b_max and b_min
    b_max = np.max([1.0, b_max])
    b_min = np.min([-1.0, b_min])

    i = y.size // 2

    if xmax is None:
        xmax = x.max()
    xmax /= scale_factor

    # Figure 1 - cross-section
    savename = f"{filename}_cross"
    fig_1, ax_1 = plt.subplots(1, 1)
    fig_1.set_size_inches(figsize)
    fig_1.set_dpi(FIG_DPI)
    fig_1.set_layout_engine("compressed")
    fig_1.suptitle("Bottom Profile - Cross-section", va="top", ha="left", x=0.01)

    ax_1.plot(
        x / scale_factor,
        b[i, :],
        rasterized=False,
        label=f"$y={y[i]/scale_factor:0.0f}$ {unit}",
    )
    _ylims = ax_1.get_ylim()
    ax_1.fill_between(
        x / scale_factor,
        _ylims[0],
        b[i, :],
        alpha=0.1,
        rasterized=False,
    )

    ax_1.axhline(color="black", linewidth=1, alpha=0.5)
    ax_1.set_xlim(0, xmax)
    ax_1.set_ylim(_ylims)
    ax_1.set_xlabel(f"$x$ [{unit}]")
    ax_1.set_ylabel("Bed Level [\si{\meter}]")
    ax_1.grid()
    ax_1.legend(loc="upper right")

    fig_1.get_layout_engine().execute(fig_1)
    save_figure(fig_1, savename)
    if not keep_open:
        plt.close(fig_1)

    # Figure 2 - map
    savename = f"{filename}_contour"
    fig_2, ax_2 = plt.subplots(1, 1)
    fig_2.set_size_inches(figsize)
    fig_2.set_dpi(FIG_DPI)
    fig_2.set_layout_engine("compressed")
    fig_2.suptitle("Bottom Profile - Contours", va="top", ha="left", x=0.01)

    cont = ax_2.contourf(
        x / scale_factor,
        y / scale_factor,
        b,
        levels=np.sort(np.linspace(0, b_min, 21)),
        cmap=cmo.tools.crop(cmo.cm.topo, b_min, b_max, 0),
        vmin=b_min,
        vmax=b_max,
    )

    cbar = fig_2.colorbar(
        cont,
        ax=ax_2,
        fraction=0.1,
        aspect=25,
        pad=0.01,
    )
    cbar.set_label("Water Depth [\si{\meter}]")
    cbar.set_ticks(np.linspace(0, b_min, 6))
    cbar.set_ticklabels(
        [f"{ticklabel:0.0f}" for ticklabel in np.linspace(0, -1.0 * b_min, 6)]
    )

    ax_2.set_xlim(0, xmax)
    ax_2.set_xlabel(f"$x$ [{unit}]")
    ax_2.set_ylabel(f"$y$ [{unit}]")

    fig_2.get_layout_engine().execute(fig_2)
    save_figure(fig_2, savename)
    if not keep_open:
        plt.close(fig_2)

    # End
    t1 = time.perf_counter_ns()
    print(f"# Finished visualising in {(t1-t0)*1e-9:0.3f} seconds")

    return (fig_1, fig_2)


if __name__ == "__main__":
    print("\nRunning inside 'bathymetry.py'")

    # Define paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    bathymetry_file = f"{PATH_TEST}/bathymetry"

    # Define function for computing 'bathymetry'
    @np.vectorize
    def function(x: Numeric, y: Numeric) -> Numeric:
        return -1.0 * np.min([x / 1000.0, 200.0])

    # Create grid and compute data
    x = np.linspace(0.0, 1e6, 7, dtype=float)
    y = np.linspace(-1e7, +1e7, 9, dtype=float)
    xx, yy = np.meshgrid(x, y)
    b = function(xx, yy)

    data = convert_to_xarray(x, y, b, savename=bathymetry_file)
    del x, y, b

    # Write data to .xyb-file
    write_bathymetry(data, filename=bathymetry_file)

    # Visualise data
    plot_bathymetry(data, filename=bathymetry_file, scale="km")
    plot_bathymetry(data, filename=bathymetry_file, half_width=True)

    print("Closing 'bathymetry.py'")
