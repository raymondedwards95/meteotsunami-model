""" Scripts to make figures for help with explaining the theory parts """

import os
import sys

import cmocean as cmo
import matplotlib.pyplot as plt
import xarray as xr

create_map: bool
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    create_map = True
except:
    print("Cartopy is not installed!")
    create_map = False

fill_map: bool
try:
    import climetlab as cml

    fill_map = True
except:
    print("Climetlab is not installed!")
    fill_map = False

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
# fmt: on


# Map
def theory_figure_map(savename: str = None, bathymetry: xr.Dataset = None) -> None:
    """Create figure of the area of interest

    Options:
        `savename`: figure name
    """
    # Arguments
    if savename is None:
        savename = f"{figure_dir}/map"

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_layout_engine("compressed")
    fig.suptitle("Area of Interest", va="top", ha="left", x=0.01)

    ax = fig.add_subplot(projection=ccrs.PlateCarree())

    if bathymetry is not None:
        savename = f"{savename}-bathymetry"

        cf = ax.pcolormesh(
            data["longitude"],
            data["latitude"],
            data,
            cmap=cmo.tools.crop(cmo.cm.topo, -200, 0, 0),
            vmin=-200,
            rasterized=True,
        )

        cb = fig.colorbar(cf, ax=ax, extend="min")
        cb.set_label("Water Depth [m]")

        ax.tick_params(labelright=False)

    coastlines = ax.coastlines("10m")
    landsurface = ax.add_feature(
        cfeature.NaturalEarthFeature(
            "physical", "land", "10m", edgecolor="black", facecolor="gray"
        )
    )
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)

    ax.set_xlim([-10, 10])
    ax.set_ylim([45, 60])
    ax.set_aspect(1)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    coastlines.set(rasterized=True)
    landsurface.set(rasterized=True)

    fig.get_layout_engine().execute(fig)
    save_figure(fig, savename)

    print(f"Saved figure {savename}")
    return


if __name__ == "__main__":
    # Settings
    show_figures = False
    current_dir = os.path.dirname(os.path.realpath(__file__))
    figure_dir = f"{current_dir}/figures"

    # Download data if necessary
    if fill_map:
        source = cml.load_source(
            "cds",
            "reanalysis-era5-single-levels",
            variable=["model_bathymetry"],
            product_type="reanalysis",
            area=[60, -15, 45, 15],
            year="2020",
            day="01",
            month="01",
            time="12:00",
            format="netcdf",
        )

        data = source.to_xarray()
        data = data["wmb"].sel(time=data["time"].min())
        data = -1.0 * data
        data = data.fillna(0.0)

        interpolate = 5
        data = data.interp(
            longitude=np.linspace(
                data["longitude"].min(),
                data["longitude"].max(),
                (data["longitude"].size - 1) * interpolate + 1,
            ),
            latitude=np.linspace(
                data["latitude"].min(),
                data["latitude"].max(),
                (data["latitude"].size - 1) * interpolate + 1,
            ),
        )

    # Parameters

    # Figures
    if create_map:
        theory_figure_map()
    if create_map and fill_map:
        theory_figure_map(bathymetry=data)

    # End
    if show_figures:
        plt.show()
