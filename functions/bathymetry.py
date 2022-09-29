""" Functions to write bathymetry data in xyz- or xyb-format

Main functions:
    convert_to_xarray
    write_bathymetry
    plot_bathymetry
"""

import os
import sys
import time
from typing import Tuple

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
        raise ValueError(
            f"'x' should be 1-dimensional, instead of {np.ndim(x)}"
        )
    if np.ndim(y) != 1:
        raise ValueError(
            f"'y' should be 1-dimensional, instead of {np.ndim(y)}"
        )
    if np.ndim(b) != 1:
        raise ValueError(
            f"'b' should be 1-dimensional, instead of {np.ndim(b)}"
        )
    raise NotImplementedError()


def _convert_to_xarray_2d(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    b: npt.ArrayLike,
) -> xr.DataArray:
    if np.ndim(x) != 1:
        raise ValueError(
            f"'x' should be 1-dimensional, instead of {np.ndim(x)}"
        )
    if np.ndim(y) != 1:
        raise ValueError(
            f"'y' should be 1-dimensional, instead of {np.ndim(y)}"
        )
    if np.ndim(b) != 2:
        raise ValueError(
            f"'b' should be 2-dimensional, instead of {np.ndim(b)}"
        )
    return xr.DataArray(b, dims=list("yx"), coords={"x": x, "y": y})


def convert_to_xarray(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    b: npt.ArrayLike,
    savename: str = None,
    close: bool = False,
) -> xr.DataArray | None:
    """ Function to convert separate numpy arrays for x, y and b to a single DataArray for writing to files

    Input:
        x:          1d array with x coordinates in km (shape = Nx)
        y:          1d array with y coordinates in km (shape = Ny)
        b:          1d or 2d array with bathymetry data in m (shape = (Ny * Nx) OR shape = (Ny, Nx))

    Options:
        savename:   saves data as an .nc file
        close:      closes data, if 'close=False', then data will be returned

    Output:
        data:       data-array containing data; if 'close=True' then data is None
    """
    print(f"# Converting bathymetry data to a data-array")
    t0 = time.perf_counter_ns()

    if b.ndim == 1:
        data = _convert_to_xarray_1d(x, y, b)
    elif b.ndim == 2:
        data = _convert_to_xarray_2d(x, y, b)
    else:
        raise ValueError(
            f"{b.ndim} dimensions is not supported for 'b' (max dim is 2"
        )

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
            encoding={"__xarray_dataarray_variable__": {
                "zlib": True,
                "complevel": 1,
                "least_significant_digit": 3
            }}
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
    """ Function to write bathymetry data to a `bathymetry.xyb` file for use with Delft3D-FM

    Input:
        data:       bathymetry and coordinate data
        filename:   .xyb file name
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

    # write
    with open(filename, "w") as file:
        file.write("")

        # loop over x
        for i in range(x.size):
            _x = x[i]

            # loop over y
            for j in range(y.size):

                # write a set of x y b on one line + remove trailing zeros
                file.write(
                    f"{_x:0.2f} {y[j]:0.2f} {b[j, i]:0.2f}\n".replace(
                        ".00", "")
                )

    t1 = time.perf_counter_ns()
    print(f"# Finished writing to '{filename}' in {(t1-t0)*1e-9:0.3f} seconds")


def plot_bathymetry(
    data: xr.DataArray,
    filename: str = None,
    xmax: Numeric = None,
    keep_open: bool = False,
) -> Tuple[plt.Figure]:
    """ Function to visualize bathymetry data

    Input:
        data:       bathymetry and coordinate data

    Options:
        filename:   name of figures
        xmax:       upper limit of x
        keep_open:  keep figures open after finishing
    """
    # Prepare
    t0 = time.perf_counter_ns()
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if filename is None:
        script_path = os.path.dirname(os.path.realpath(__file__))
        filename = f"{script_path}/tests/fig_bathymetry"
    if filename.endswith(".jpg"):
        filename.replace(".jpg", "")
    print(f"# Visualizing bathymetry in '{filename}'")

    # Extract data
    x = data.x.values
    y = data.y.values
    b = data.values

    b_max = np.ceil(b.max())
    b_min = np.floor(b.min())

    # fix b_max and b_min
    b_max = np.max([1., b_max])
    b_min = np.min([-1., b_min])

    i = y.size // 2

    if xmax is None:
        xmax = x.max()
    xmax /= 1000

    # Figure 1 - cross-section
    savename = f"{filename}_cross"
    fig_1, ax_1 = plt.subplots(1, 1)
    fig_1.set_size_inches(FIGSIZE_NORMAL)
    fig_1.set_dpi(FIG_DPI)
    fig_1.set_layout_engine("compressed")
    ax_1.plot(x / 1000., b[i, :])
    _ylims = ax_1.get_ylim()
    ax_1.fill_between(x / 1000., _ylims[0], b[i, :], alpha=0.1)
    ax_1.set_title(f"Bottom Profile Cross-Section at $y={y[i]}$")
    ax_1.set_xlim(0, xmax)
    ax_1.set_ylim(_ylims)
    ax_1.set_xlabel("$x$ [km]")
    ax_1.set_ylabel("Bed Level [m]")
    ax_1.grid()
    fig_1.get_layout_engine().execute(fig_1)
    fig_1.savefig(savename, bbox_inches="tight",
                  dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)
    if not keep_open:
        plt.close(fig_1)
    print(f"# Saved figure '{savename}'")

    # Figure 2 - map
    savename = f"{filename}_contour"
    fig_2, ax_2 = plt.subplots(1, 1)
    fig_2.set_size_inches(FIGSIZE_NORMAL)
    fig_2.set_dpi(FIG_DPI)
    fig_2.set_layout_engine("compressed")
    div = make_axes_locatable(ax_2)
    cax = div.append_axes("right", "5%", "5%")
    cont = ax_2.contourf(
        x / 1000.,
        y / 1000.,
        b,
        levels=np.sort(np.linspace(0, b_min, 21)),
        cmap=cmo.tools.crop(cmo.cm.topo, b_min, b_max, 0),
        vmin=b_min,
        vmax=b_max,
    )
    cbar = fig_2.colorbar(cont, cax=cax)
    cbar.set_label("Water Depth [m]")
    cbar.set_ticks(np.linspace(0, b_min, 6))
    cbar.set_ticklabels(
        [f"{ticklabel:0.0f}" for ticklabel in np.linspace(0, -1. * b_min, 6)])
    ax_2.set_title(f"Bottom Profile")
    ax_2.set_xlim(0, xmax)
    ax_2.set_xlabel("$x$ [km]")
    ax_2.set_ylabel("$y$ [km]")
    fig_2.get_layout_engine().execute(fig_2)
    fig_2.savefig(savename, bbox_inches="tight",
                  dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)
    if not keep_open:
        plt.close(fig_2)
    print(f"# Saved figure '{savename}'")

    # End
    t1 = time.perf_counter_ns()
    print(f"# Finished visualising in {(t1-t0)*1e-9:0.3f} seconds")

    return (fig_1, fig_2)


if __name__ == "__main__":
    # Define paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    bathymetry_dir = f"{script_dir}/tests/bath"
    bathymetry_file = f"{bathymetry_dir}/bathymetry"
    os.makedirs(bathymetry_dir, exist_ok=True)

    # Define function for computing 'bathymetry'
    @np.vectorize
    def function(x: Numeric, y: Numeric) -> Numeric:
        return -1. * np.min([x / 1000., 200.])

    # Create grid and compute data
    x = np.linspace(0., 1e6, 5, dtype=float)
    y = np.linspace(-1e7, +1e7, 5, dtype=float)
    xx, yy = np.meshgrid(x, y)
    b = function(xx, yy)

    data = convert_to_xarray(x, y, b, savename=bathymetry_file)
    del x, y, b

    # Write data to .xyb-file
    write_bathymetry(data, filename=bathymetry_file)

    # Visualise data
    plot_bathymetry(data, filename=bathymetry_file)
