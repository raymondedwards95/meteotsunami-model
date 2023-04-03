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
import functions.utilities as fu
# fmt: on


# Data
LOCATIONS = {
    "Vlissingen": {
        "lon": 3.566667,
        "lat": 51.45,
    },
    "Hoek van Holland": {
        "lon": 4.128611,
        "lat": 51.981111,
    },
    "Den Helder": {
        "lon": 4.75,
        "lat": 52.933333,
    },
}


# Map
def theory_figure_map(
    savename: str = None,
    bathymetry: xr.Dataset = None,
    locations: bool = False,
    zoom: bool = False,
) -> None:
    """Create figure of the area of interest

    Options:
        `savename`:     figure name
        `bathymetry`:   plot data from given Dataset
        `locations`:    plot points for specific locations if true
        `zoom`:         focus on small part of full area
    """
    max_depth = -200.0

    # Arguments
    if savename is None:
        savename = f"{figure_dir}/map"

    if zoom:
        savename = f"{savename}-zoom"
        max_depth = -50.0

    # Figure
    fig = plt.figure()
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_layout_engine("compressed")
    # fig.suptitle("Area of Interest")

    ax = fig.add_subplot(projection=ccrs.PlateCarree())

    if bathymetry is not None:
        savename = f"{savename}-bathymetry"

        cf = ax.pcolormesh(
            data["longitude"],
            data["latitude"],
            data,
            cmap=cmo.tools.crop(cmo.cm.topo, max_depth, 0, 0),
            vmin=max_depth,
            rasterized=True,
        )

        cb = fig.colorbar(cf, ax=ax, extend="min")
        cb.set_label("Bed Level [\\si{\\meter}]")

        ax.tick_params(labelright=False)

    if locations:
        for i, loc in enumerate(["Vlissingen", "Hoek van Holland", "Den Helder"]):
            ax.plot(
                [LOCATIONS[loc]["lon"]],
                [LOCATIONS[loc]["lat"]],
                "*",
                label=f"{loc}",
                color=f"C{i}",
                markersize=10,
                markeredgecolor="black",
                markeredgewidth=0.8,
            )
            ax.annotate(
                text=f"{loc}",
                xy=[LOCATIONS[loc]["lon"], LOCATIONS[loc]["lat"]],
                xytext=[0.98, 0.03 + i * 0.08],
                textcoords="axes fraction",
                color=f"C{i}",
                ha="right",
                va="bottom",
                arrowprops=dict(
                    arrowstyle="->",
                    connectionstyle="angle,angleA=180,angleB=90",
                    shrinkA=0,
                    shrinkB=5,
                ),
                bbox=dict(
                    facecolor="white",
                    edgecolor="white",
                ),
            )

    coastlines = ax.coastlines("10m")
    landsurface = ax.add_feature(
        cfeature.NaturalEarthFeature(
            "physical",
            "land",
            "10m",
            edgecolor="black",
            facecolor="silver",
        )
    )
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)

    if zoom:
        ax.set_xlim([1.0, 6.0])
        ax.set_ylim([50.5, 54.5])
    else:
        ax.set_xlim([-10.0, 10.0])
        ax.set_ylim([45.0, 60.0])
    ax.set_aspect(1)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    coastlines.set(rasterized=True)
    landsurface.set(rasterized=True)

    fig.get_layout_engine().execute(fig)
    save_figure(fig, savename)

    print(f"Saved figure {savename}")
    return


# Cross-section bathymetry
def theory_figure_cross(
    bathymetry: xr.DataArray,
    savename: str = None,
) -> None:
    """Create figure of cross sections of bathymetry

    Input:
        `bathymetry`:   data with bathymetry

    Options:
        `savename`:     figure name
    """
    # Arguments
    if savename is None:
        savename = f"{figure_dir}/map-cross-section"

    if isinstance(bathymetry, xr.Dataset):
        raise TypeError(
            f"Parameter `bathymetry` should be xr.DataArray instead of `{type(bathymetry)}`"
        )

    dist_max = 2e5

    # Figure
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_layout_engine("compressed")
    # fig.suptitle("Cross-sections of bottom topography")

    for i, loc in enumerate(["Vlissingen", "Hoek van Holland", "Den Helder"]):
        # Select longitude
        lon_max = LOCATIONS[loc]["lon"]
        lon_min = LOCATIONS[loc]["lon"] - 5
        lons = bathymetry["longitude"].sel(longitude=slice(lon_min, lon_max)).values
        # lons = np.linspace(lon_min, lon_max, 500)
        lons = xr.DataArray(lons, dims="distance")

        # Select latitude
        lat_max = LOCATIONS[loc]["lat"] + 5
        lat_min = LOCATIONS[loc]["lat"]
        lats = bathymetry["latitude"].sel(latitude=slice(lat_min, lat_max)).values[::-1]
        # lats = np.linspace(lat_min, lat_max, 500)[::-1]
        lats = xr.DataArray(lats, dims="distance")

        # Interpolate over section
        cross = bathymetry.interp(latitude=lats, longitude=lons)

        # Compute distance from reference
        cross["distance"] = fu.haversine_distance(
            LOCATIONS[loc]["lon"],
            LOCATIONS[loc]["lat"],
            lons.values,
            lats.values,
        )

        # Selecting data by distance
        cross = cross.where(cross["distance"] <= dist_max)

        # Plot
        ax.plot(
            cross["distance"] / 1e3,
            cross,
            label=loc,
            rasterized=False,
        )
        # ax.fill_between(
        #     cross["distance"] / 1e3,
        #     cross,
        #     alpha=0.1,
        #     rasterized=False,
        # )

    ax.legend()
    ax.grid()
    ax.set_xlabel("Distance from Location [\\si{\\kilo\\meter}]")
    ax.set_ylabel("Bed Level [\\si{\\meter}]")
    ax.axhline(color="black", linewidth=1, alpha=0.5)
    ax.set_xlim(0.0, dist_max / 1e3)

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

        data = source.to_xarray(chunks="auto")
        data = data["wmb"].sel(time=data["time"].min()).drop_vars("time")
        data = data.astype(float)
        data = -1.0 * data
        data = data.fillna(0.0)

        interpolate = 5
        if interpolate > 1:
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
        theory_figure_map(zoom=True)
    if create_map and fill_map:
        theory_figure_map(bathymetry=data, locations=True)
        theory_figure_map(bathymetry=data, locations=True, zoom=True)
    if fill_map:
        theory_figure_cross(bathymetry=data)

    # End
    if show_figures:
        plt.show()
