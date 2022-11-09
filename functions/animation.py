""" Functions for visualising model outputs

Main functions:
    animation_contour
    animation_contour_uv
    animation_alongshore
    animation_crossshore
"""

import os
import sys
import time
import warnings

import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
# fmt: on


def animation_contour(
    data: xr.Dataset,
    savedir: str,
    xlims: Numeric = None,
    _test_i_max: Integer = None,
    close: bool = True,
    fps: Numeric = 20,
    static: bool = False,
) -> None:
    """Creates an animation of the top-down view of the water level and surface air pressure data

    Input:
        `data`:     Dataset that contains all variables
        `savedir`:  Folder to write the animation to

    Options:
        `xlims`:    x-limits for the figure
        `close`:    close figure after finishing
    """
    # Check inputs
    if not isinstance(data, xr.Dataset):
        raise TypeError(f"data is not a Dataset, but it is {type(data)}")

    os.makedirs(savedir, exist_ok=True)
    savename = f"{savedir}_contours.{ANIM_EXT}"
    savename_static = f"{savedir}__static_contours"

    # Shortcuts
    x = data["x"]
    y = data["y"]
    t = data["t"]
    wl = data["wl"]
    p = data["p"]

    x = x - x.min()
    wl_max = float(np.ceil(np.max([np.abs([wl.max(), wl.min()])])))
    p_max = float(np.ceil(np.max([np.abs([p.max(), p.min()])])))

    limits = [(-1.0 * wl_max, wl_max), (-1.0 * p_max, p_max)]
    cmaps = [cmo.cm.balance, cmo.cm.curl]

    wl_levels = np.linspace(
        np.floor(wl.min() * 10.0) / 10.0,
        np.ceil(wl.max() * 10.0) / 10.0,
        101,
    )
    wl_ticks = np.linspace(
        np.floor(wl.min() * 10.0) / 10.0,
        np.ceil(wl.max() * 10.0) / 10.0,
        5,
    )
    p_levels = np.linspace(
        np.floor(p.min() * 10.0) / 10.0,
        np.ceil(p.max() * 10.0) / 10.0,
        101,
    )
    p_ticks = np.linspace(
        np.floor(p.min() * 10.0) / 10.0,
        np.ceil(p.max() * 10.0) / 10.0,
        5,
    )

    if xlims is None:
        xlims = [y.min() / 1000.0 / 10.0, y.max() / 1000.0]

    # Figure options
    fig, ax = plt.subplots(2, 1, sharey=True)
    fig.set_dpi(200)
    fig.set_size_inches(1440 / fig.get_dpi(), 720 / fig.get_dpi())
    fig.set_tight_layout(True)

    div = np.array([make_axes_locatable(ax[i]) for i in range(2)])
    cax = np.array([div[i].append_axes("right", "5%", "5%") for i in range(2)])
    cbar = np.array([None, None])

    # Initial data
    plotdata = np.zeros(2, dtype=object)
    plottext = np.zeros(1, dtype=object)

    def set_plotdata(i=0):
        plotdata[0] = ax[0].contourf(
            y / 1000,
            x / 1000,
            wl.isel(t=i).T,
            vmin=-1.0 * wl_max,
            vmax=wl_max,
            cmap=cmaps[0],
            levels=wl_levels,
        )
        plotdata[1] = ax[1].contourf(
            y / 1000,
            x / 1000,
            p.isel(t=i).T,
            vmin=-1.0 * p_max,
            vmax=p_max,
            cmap=cmaps[1],
            levels=p_levels,
        )

    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )

    set_plotdata()
    set_plottext()

    # Subplot layout
    def initfig():
        for i in range(2):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_ylabel("$x$ [km]")
            ax[i].set_ylim([0, 200])
            ax[i].set_xlim(xlims)

            # colorbar
            plotdata[i].set_clim(limits[i])
            cax[i].cla()
            cbar[i] = fig.colorbar(plotdata[i], cax=cax[i])

        cbar[0].set_ticks(wl_ticks)
        cbar[1].set_ticks(p_ticks)

        ax[0].set_xlabel("$y$ [km]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    initfig()

    # Update data
    num_frames = t.size
    if _test_i_max is not None and isinstance(_test_i_max, int):
        warnings.warn("Parameter _test_i_max is set for testing animation function!")
        num_frames = _test_i_max

    def update(i):
        # progress
        if not (num_frames - i - 1) % (num_frames // 5):
            t0_interval = time.perf_counter_ns()
            print(
                f"# Frame {i+1:4.0f} of {num_frames:0.0f} ({(i+1)/num_frames*100:0.1f}%) ({(i+1) / ((t0_interval-t0) * 1e-9):0.1f} fps)"
            )

        # remove data in contour plots
        for _temp in plotdata[0].collections:
            _temp.remove()
        for _temp in plotdata[1].collections:
            _temp.remove()

        # new data
        set_plotdata(i)
        set_plottext(i)

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    # Animation
    frames = (np.arange(num_frames)).astype(int)
    anim = FuncAnimation(
        fig,
        update,
        init_func=initfig,
        frames=frames,
        interval=1000 / fps,
    )
    print(f"# Creating animation '{savename}'")
    t0 = time.perf_counter_ns()
    if static: plt.savefig(savename_static)
    anim.save(savename, extra_args=FFMPEG_ARGS)
    t1 = time.perf_counter_ns()
    print(
        f"# Finished contour-animation in {(t1-t0) * 1e-9:0.1f} seconds (average {num_frames / ((t1-t0) * 1e-9):0.1f} frames per second)"
    )
    print(f"# Saved animation as '{savename}'")

    if close:
        plt.close(fig)
        print(f"# Figure is closed")
    else:
        print(f"# Figure is not closed!")

    return


def animation_contour_uv(
    data: xr.Dataset,
    savedir: str,
    xlims: Numeric = None,
    _test_i_max: Integer = None,
    close: bool = True,
    fps: Numeric = 20,
    static: bool = False,
) -> None:
    """Creates an animation of the top-down view of the water velocity data

    Input:
        `data`:     Dataset that contains all variables
        `savedir`:  Folder to write the animation to

    Options:
        `xlims`:    x-limits for the figure
        `close`:    close figure after finishing
    """
    # Check inputs
    if not isinstance(data, xr.Dataset):
        raise TypeError(f"data is not a Dataset, but it is {type(data)}")

    os.makedirs(savedir, exist_ok=True)
    savename = f"{savedir}_uv_cont.{ANIM_EXT}"
    savename_static = f"{savedir}__static_uv_cont"

    # Shortcuts
    x = data["x"]
    y = data["y"]
    t = data["t"]
    u = data["u"]
    v = data["v"]

    x = x - x.min()
    uv_max = float(
        np.ceil(10.0 * np.max([np.abs([u.max(), u.min(), v.max(), v.min()])])) / 10.0
    )

    limits = (-1.0 * uv_max, uv_max)
    cmaps = [cmo.cm.delta, cmo.cm.delta]

    if xlims is None:
        xlims = [y.min() / 1000.0 / 10.0, y.max() / 1000.0]

    uv_levels = np.linspace(-1.0 * uv_max, uv_max, 101)
    uv_ticks = np.linspace(-1.0 * uv_max, uv_max, 5)

    # Figure options
    fig, ax = plt.subplots(2, 1, sharey=True)
    fig.set_dpi(200)
    fig.set_size_inches(1440 / fig.get_dpi(), 720 / fig.get_dpi())
    fig.set_tight_layout(True)

    div = np.array([make_axes_locatable(ax[i]) for i in range(2)])
    cax = np.array([div[i].append_axes("right", "5%", "5%") for i in range(2)])
    cbar = np.array([None, None])

    # Initial data
    plotdata = np.zeros(2, dtype=object)
    plottext = np.zeros(1, dtype=object)

    def set_plotdata(i=0):
        plotdata[0] = ax[0].contourf(
            y / 1000,
            x / 1000,
            u.isel(t=i).T,
            vmin=-1.0 * uv_max,
            vmax=uv_max,
            cmap=cmaps[0],
            levels=uv_levels,
        )
        plotdata[1] = ax[1].contourf(
            y / 1000,
            x / 1000,
            v.isel(t=i).T,
            vmin=-1.0 * uv_max,
            vmax=uv_max,
            cmap=cmaps[1],
            levels=uv_levels,
        )

    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )

    set_plotdata()
    set_plottext()

    # Subplot layout
    def initfig():
        for i in range(2):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_ylabel("$x$ [km]")
            ax[i].set_ylim([0, 200])
            ax[i].set_xlim(xlims)

            # colorbar
            plotdata[i].set_clim(limits)
            cax[i].cla()
            cbar[i] = fig.colorbar(plotdata[i], cax=cax[i])

        cbar[0].set_ticks(uv_ticks)
        cbar[1].set_ticks(uv_ticks)

        ax[0].set_xlabel("$y$ [km]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    initfig()

    # Update data
    num_frames = t.size
    if _test_i_max is not None and isinstance(_test_i_max, int):
        warnings.warn("Parameter _test_i_max is set for testing animation function!")
        num_frames = _test_i_max

    def update(i):
        # progress
        if not (num_frames - i - 1) % (num_frames // 5):
            t0_interval = time.perf_counter_ns()
            print(
                f"# Frame {i+1:4.0f} of {num_frames:0.0f} ({(i+1)/num_frames*100:0.1f}%) ({(i+1) / ((t0_interval-t0) * 1e-9):0.1f} fps)"
            )

        # remove data in contour plots
        for _temp in plotdata[0].collections:
            _temp.remove()
        for _temp in plotdata[1].collections:
            _temp.remove()

        # new data
        set_plotdata(i)
        set_plottext(i)

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    # Animation
    frames = (np.arange(num_frames)).astype(int)
    anim = FuncAnimation(
        fig,
        update,
        init_func=initfig,
        frames=frames,
        interval=1000 / fps,
    )
    print(f"# Creating animation '{savename}'")
    t0 = time.perf_counter_ns()
    if static: plt.savefig(savename_static)
    anim.save(savename, extra_args=FFMPEG_ARGS)
    t1 = time.perf_counter_ns()
    print(
        f"# Finished uv-contour-animation in {(t1-t0) * 1e-9:0.1f} seconds (average {num_frames / ((t1-t0) * 1e-9):0.1f} frames per second)"
    )
    print(f"# Saved animation as '{savename}'")

    if close:
        plt.close(fig)
        print(f"# Figure is closed")
    else:
        print(f"# Figure is not closed!")

    return


def animation_alongshore(
    data: xr.Dataset,
    savedir: str,
    xlims: Numeric = None,
    _test_i_max: Integer = None,
    close: bool = True,
    fps: Numeric = 20,
    static: bool = False,
) -> None:
    """Creates an animation of an alongshore cross-section of the water level and surface air pressure data

    Input:
        `data`:     Dataset that contains all variables
        `savedir`:  Folder to write the animation to

    Options:
        `xlims`:    x-limits for the figure
        `close`:    close figure after finishing
    """
    # Check inputs
    if not isinstance(data, xr.Dataset):
        raise TypeError(f"data is not a Dataset, but it is {type(data)}")

    os.makedirs(savedir, exist_ok=True)
    savename = f"{savedir}_alongshore.{ANIM_EXT}"
    savename_static = f"{savedir}__static_alongshore"

    # Shortcuts
    x = data["x"]
    y = data["y"]
    t = data["t"]
    wl = data["wl"]
    p = data["p"]

    x = x - x.min()
    wl_max = float(np.max([np.abs([wl.max(), wl.min()])]))
    p_max = float(p.max())
    p_min = float(p.min())

    if xlims is None:
        xlims = [y.min() / 1000.0 / 10.0, y.max() / 1000.0]

    # Figure options
    fig, ax = plt.subplots(2, 1, sharex=True)
    fig.set_dpi(200)
    fig.set_size_inches(1440 / fig.get_dpi(), 720 / fig.get_dpi())
    fig.set_tight_layout(True)

    # Initial data
    plotdata = np.zeros(2, dtype=object)
    plottext = np.zeros(1, dtype=object)

    def set_plotdata(i=0):
        plotdata[0] = ax[0].plot(y / 1000, wl.isel(t=i).interp(x=10000), color="C0")[0]
        plotdata[1] = ax[1].plot(y / 1000, p.isel(t=i).interp(x=10000), color="C1")[0]

    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )

    set_plotdata()
    set_plottext()

    # Subplot layout
    def initfig():
        for i in range(2):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_xlim(xlims)

        ax[0].set_ylim([-1.0 * wl_max, wl_max])
        ax[1].set_ylim([p_min, p_max])
        ax[1].set_xlabel("$y$ [km]")
        ax[0].set_ylabel("Sea Surface Elevation [m]")
        ax[1].set_ylabel("Surface Air Pressure [Pa]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    initfig()

    # Update data
    num_frames = t.size
    if _test_i_max is not None and isinstance(_test_i_max, int):
        warnings.warn("Parameter _test_i_max is set for testing animation function!")
        num_frames = _test_i_max

    def update(i):
        # progress
        if not (num_frames - i - 1) % (num_frames // 5):
            t0_interval = time.perf_counter_ns()
            print(
                f"# Frame {i+1:4.0f} of {num_frames:0.0f} ({(i+1)/num_frames*100:0.1f}%) ({(i+1) / ((t0_interval-t0) * 1e-9):0.1f} fps)"
            )

        # new data
        plotdata[0].set_ydata(wl.isel(t=i).interp(x=10000))
        plotdata[1].set_ydata(p.isel(t=i).interp(x=10000))
        set_plottext(i)

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    # Animation
    frames = (np.arange(num_frames)).astype(int)
    anim = FuncAnimation(
        fig,
        update,
        init_func=initfig,
        frames=frames,
        interval=1000 / fps,
    )
    print(f"# Creating animation '{savename}'")
    t0 = time.perf_counter_ns()
    if static: plt.savefig(savename_static)
    anim.save(savename, extra_args=FFMPEG_ARGS)
    t1 = time.perf_counter_ns()
    print(
        f"# Finished alongshore-animation in {(t1-t0) * 1e-9:0.1f} seconds (average {num_frames / ((t1-t0) * 1e-9):0.1f} frames per second)"
    )
    print(f"# Saved animation as '{savename}'")

    if close:
        plt.close(fig)
        print(f"# Figure is closed")
    else:
        print(f"# Figure is not closed!")

    return


def animation_crossshore(
    data: xr.Dataset,
    savedir: str,
    _test_i_max: Integer = None,
    close: bool = True,
    fps: Numeric = 20,
    static: bool = False,
) -> None:
    """Creates an animation of crossshore cross-sections of the water level and surface air pressure data

    Input:
        `data`:     Dataset that contains all variables
        `savedir`:  Folder to write the animation to

    Options:
        `close`:    close figure after finishing
    """
    # Check inputs
    if not isinstance(data, xr.Dataset):
        raise TypeError(f"data is not a Dataset, but it is {type(data)}")

    os.makedirs(savedir, exist_ok=True)
    savename = f"{savedir}_crossshore.{ANIM_EXT}"
    savename_static = f"{savedir}__static_crossshore"

    # Shortcuts
    x = data["x"]
    y = data["y"]
    t = data["t"]
    wl = data["wl"]

    x = x - x.min()
    wl_max = float(np.max([np.abs([wl.max(), wl.min()])]))
    wl_min = -1.0 * wl_max
    y_max = float(y.max())

    # Figure options
    slices = 5
    fig, ax = plt.subplots(slices, 1, sharex=True)
    fig.set_dpi(200)
    fig.set_size_inches(1440 / fig.get_dpi(), 720 / fig.get_dpi())
    fig.set_tight_layout(True)

    # Initial data
    plotdata = np.zeros(slices, dtype=object)
    plottext = np.zeros(1, dtype=object)
    _y = np.arange(slices) * y_max / slices

    def set_plotdata(i=0):
        for j in range(slices):
            plotdata[j] = ax[j].plot(
                x / 1000,
                wl.isel(t=i).interp(y=_y[j]),
                color="C0",
                label=f"$y={_y[j]/1000:0.0f}$ km",
            )[0]

    def set_plottext(i=0):
        plottext[0] = ax[0].set_title(
            f"$t={t.isel(t=i).values.tolist()/1e9/3600:0.1f}$ hours since start"
        )

    set_plotdata()
    set_plottext()

    # Subplot layout
    def initfig():
        for i in range(slices):
            ax[i].axhline(color="black", linewidth=1)
            ax[i].axvline(color="black", linewidth=1)
            ax[i].set_xlim([0, 200])
            ax[i].set_ylabel("SSE [m]")
            ax[i].set_ylim([wl_min, wl_max])
            ax[i].legend()

        ax[-1].set_xlabel("$x$ [km]")
        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    initfig()

    # Update data
    num_frames = t.size
    if _test_i_max is not None and isinstance(_test_i_max, int):
        warnings.warn("Parameter _test_i_max is set for testing animation function!")
        num_frames = _test_i_max

    def update(i):
        # progress
        if not (num_frames - i - 1) % (num_frames // 5):
            t0_interval = time.perf_counter_ns()
            print(
                f"# Frame {i+1:4.0f} of {num_frames:0.0f} ({(i+1)/num_frames*100:0.1f}%) ({(i+1) / ((t0_interval-t0) * 1e-9):0.1f} fps)"
            )

        # new data
        for j in range(slices):
            plotdata[j].set_ydata(wl.isel(t=i).interp(y=_y[j]))
        set_plottext(i)

        return tuple(plotdata.flatten()) + tuple(plottext.flatten())

    # Animation
    frames = (np.arange(num_frames)).astype(int)
    anim = FuncAnimation(
        fig,
        update,
        init_func=initfig,
        frames=frames,
        interval=1000 / fps,
    )
    print(f"# Creating animation '{savename}'")
    t0 = time.perf_counter_ns()
    if static: plt.savefig(savename_static)
    anim.save(savename, extra_args=FFMPEG_ARGS)
    t1 = time.perf_counter_ns()
    print(
        f"# Finished crossshore-animation in {(t1-t0) * 1e-9:0.1f} seconds (average {num_frames / ((t1-t0) * 1e-9):0.1f} frames per second)"
    )
    print(f"# Saved animation as '{savename}'")

    if close:
        plt.close(fig)
        print(f"# Figure is closed")
    else:
        print(f"# Figure is not closed!")

    return


if __name__ == "__main__":
    # Define paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    anim_dir = f"{script_dir}/tests/anim/anim"
    os.makedirs(anim_dir, exist_ok=True)

    # Read data
    data = xr.open_dataset(
        f"{script_dir}/../reproduction-an-2012/output/data_repr_17.nc",
        # chunks={"y": -1, "x": -1, "t": "auto"},
    )
    print(f"Grid size: {data.sizes}")
    # print(f"Chunksize: {data.chunksizes}")

    # Make animations
    animation_contour(data, savedir=anim_dir, _test_i_max=25, static=True,)
    animation_contour_uv(data, savedir=anim_dir, _test_i_max=25, static=True,)

    animation_alongshore(data, savedir=anim_dir, _test_i_max=25, static=True,)
    animation_crossshore(data, savedir=anim_dir, _test_i_max=15, fps=10, static=True,)

    # End
    plt.close("all")
    data.close()
