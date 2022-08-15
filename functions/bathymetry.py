""" Functions to write bathymetry data in xyz- or xyb-format """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


def _convert_to_xarray_1d(x: npt.ArrayLike, y: npt.ArrayLike, b: npt.ArrayLike) -> xr.DataArray:
    if np.ndim(x) != 1:
        raise ValueError(f"'x' should be 1-dimensional, instead of {np.ndim(x)}")
    if np.ndim(y) != 1:
        raise ValueError(f"'y' should be 1-dimensional, instead of {np.ndim(x)}")
    if np.ndim(b) != 1:
        raise ValueError(f"'b' should be 1-dimensional, instead of {np.ndim(x)}")
    raise NotImplementedError()


def _convert_to_xarray_2d(x: npt.ArrayLike, y: npt.ArrayLike, b: npt.ArrayLike) -> xr.DataArray:
    if np.ndim(x) != 1:
        raise ValueError(f"'x' should be 1-dimensional, instead of {np.ndim(x)}")
    if np.ndim(y) != 1:
        raise ValueError(f"'y' should be 1-dimensional, instead of {np.ndim(x)}")
    if np.ndim(b) != 2:
        raise ValueError(f"'b' should be 2-dimensional, instead of {np.ndim(x)}")
    return xr.DataArray(b, dims=list("yx"), coords={"x": x, "y": y})


def convert_to_xarray(x: npt.ArrayLike, y: npt.ArrayLike, b: npt.ArrayLike, savename: str=None, close: bool=False) -> xr.DataArray | None:
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
    if b.ndim == 1:
        data = _convert_to_xarray_1d(x, y, b)
    elif b.ndim == 2:
        data = _convert_to_xarray_2d(x, y, b)
    else:
        raise ValueError(f"{b.ndim} dimensions is not supported for 'b' (max dim is 2")

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
        print(f"Saved data-array as {savename}")

    if close:
        return
    return data


def write_bathymetry(data: xr.DataArray, filename) -> None:
    """ Function to write bathymetry data to a `bathymetry.xyb` file for use with Delft3D-FM

    Input:
        data:       bathymetry and coordinate data
        filename:   .xyb file name
    """
    ## prepare
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if not filename.endswith(".xyb"):
        filename += ".xyb"
    print(f"Writing bathymetry data to '{filename}'")

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


def plot_bathymetry(data: xr.DataArray, filename: str=None, xmax: Numeric=None, keep_open: bool=False) -> None:
    """ Function to visualize bathymetry data

    Input:
        data:       bathymetry and coordinate data

    Options:
        filename:   name of figures
        xmax:       upper limit of x
        keep_open:  keep figures open after finishing
    """
    ## Prepare
    assert type(data) == xr.DataArray, "Input is not a DataArray"

    if filename is None:
        filename = os.path.dirname(os.path.realpath(__file__)) + "/tests/fig_bathymetry"
    if filename.endswith(".jpg"):
        filename.replace(".jpg", "")
    print(f"Visualizing bathymetry in '{filename}'")

    ## Extract data
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

    ## Figure 1 - cross-section
    savename = f"{filename}_cross"
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    ax.plot(x / 1000., b[i, :])
    _ylims = ax.get_ylim()
    ax.fill_between(x / 1000., _ylims[0], b[i, :], alpha=0.1)
    ax.set_title(f"Bottom Profile Cross-Section at $y={y[i]}$")
    ax.set_xlim(0, xmax)
    ax.set_ylim(_ylims)
    ax.set_xlabel("$x$ [km]")
    ax.set_ylabel("Bed Level [m]")
    ax.grid()
    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    if not keep_open:
        plt.close(fig)
    print(f"Saved figure '{savename}'")

    ## Figure 2 - map
    savename = f"{filename}_contour"
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", "5%", "5%")
    cont = ax.contourf(
        x / 1000.,
        y / 1000.,
        b,
        levels=np.sort(np.linspace(0, b_min, 21)),
        cmap=cmo.tools.crop(cmo.cm.topo, b_min, b_max, 0),
        vmin=b_min,
        vmax=b_max,
    )
    cbar = fig.colorbar(cont, cax=cax)
    cbar.set_label("Water Depth [m]")
    cbar.set_ticks(np.linspace(0, b_min, 6))
    cbar.set_ticklabels([f"{ticklabel:0.0f}" for ticklabel in np.linspace(0, -1. * b_min, 6)])
    ax.set_title(f"Bottom Profile")
    ax.set_xlim(0, xmax)
    ax.set_xlabel("$x$ [km]")
    ax.set_ylabel("$y$ [km]")
    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs={"optimize": True, "compress_level": 9})
    if not keep_open:
        plt.close(fig)
    print(f"Saved figure '{savename}'")

    return


if __name__ == "__main__":
    import time

    # Paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    bathymetry_dir = f"{script_dir}/tests/bathymetry"

    # function for computing 'bathymetry'
    @np.vectorize
    def function(x: Numeric, y: Numeric) -> Numeric:
        return -1. * np.min([x / 1000., 200.])

    t0 = time.perf_counter()

    # grid and data
    x = np.linspace(0., 1e6, 5, dtype=float)
    y = np.linspace(-1e7, +1e7, 5, dtype=float)
    xx, yy = np.meshgrid(x, y)
    b = function(xx, yy)

    data = convert_to_xarray(x, y, b, savename=bathymetry_dir)
    del x, y, b

    # write data to .xyb-file
    write_bathymetry(data, filename=bathymetry_dir)

    # evaluate performance
    t1 = time.perf_counter()
    print(f"Created test bathymetry data in {t1 - t0 :0.2f} seconds.")

    plot_bathymetry(data, filename=bathymetry_dir)

    # evaluate performance
    t2 = time.perf_counter()
    print(f"Visualized test bathymetry data in {t2 - t1 :0.2f} seconds.")

    #plt.show()
