""" Functions to create observation points and cross sections for D3D-FM """

import os
import sys

import matplotlib.axes as mpa
import matplotlib.pyplot as plt
import numpy.typing as npt

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *


class ObservationPoint():
    """ Single point for observation

    Attributes:
        x:      x-coordinate
        y:      y-coordinate
        name:   (unique) name
    """

    def __init__(self, x: Numeric, y: Numeric, name: str = "") -> None:
        """ Create an observation point

        Input:
            x:      x-coordinate
            y:      y-coordinate

        Options:
            name:   (unique) name
        """
        self.x = x
        self.y = y
        self.name = name

    def __str__(self):
        return f"Observation point '{self.name}'' at x={self.x:.0e} and y={self.y:.0e}".replace("e+00", "")

    def plot(self, ax: mpa.Axes) -> mpa.Axes:
        """ Plot point on an existing Axes object

        Input:
            ax:     existing Axes object
        """
        ax.plot(self.x, self.y, "o", markersize=5, label=self.name)
        return ax


class ObservationCrossSection():
    """ Lines for observation consisting of two or more points

    Attributes:
        x:      x-coordinates of connecting points
        y:      y-coordinates of connecting points
        name:   (unique) name
        n:      number of points
    """

    def __init__(self, x: npt.ArrayLike, y: npt.ArrayLike, name: str = "") -> None:
        """ Create an observation point

        Input:
            x:      array of x-coordinates
            y:      array of y-coordinates

        Options:
            name:   (unique) name
        """
        self.x = list(x)
        self.y = list(y)
        self.n = len(x)
        self.name = name

        assert len(x) > 1
        assert len(x) == len(y)

    def __str__(self):
        return f"Observation cross-section '{self.name}' between points p_start=({self.x[0]:.0e}, {self.y[0]:.0e}) and p_end=({self.x[-1]:.0e}, {self.y[-1]:.0e})".replace("e+00", "")

    def plot(self, ax: mpa.Axes) -> mpa.Axes:
        """ Plot lines on an existing Axes object

        Input:
            ax:     existing Axes object
        """
        ax.plot(self.x, self.y, "-o", markersize=5, label=self.name)
        return ax


def write_observations(data: npt.ArrayLike, filename: str) -> None:
    """ Write observation points and observation cross sections to files

    Input:
        data:       list of observation points and observation cross sections
        filename:   name of file to write (extensions are added automatically)
    """
    print("\nWriting observation points and cross-sections")

    # Path and filenames
    if filename.endswith(".xyn"):
        filename.replace(".xyn", "")
    if filename.endswith("_crs.pli"):
        filename.replace("_crs.pli", "")

    filename_points = filename + ".xyn"
    filename_sections = filename + "_crs.pli"

    # Sort points and cross sections
    data = list(data)

    points = []
    sections = []

    for element in data:
        if isinstance(element, ObservationPoint):
            points.append(element)
        elif isinstance(element, ObservationCrossSection):
            sections.append(element)
        else:
            print(f"'{element}' is not an ObservationPoint or ObservationCrossSection")

    # Process points
    if len(points) > 0:
        with open(filename_points, "w") as file:
            for element in points:
                print(f"* {element}")
                _name = element.name.replace("\'", "")
                file.write(f"{element.x} \t{element.y} \t'{_name}'\n")

    print(f"Finished writing points to '{filename_points}'")

    # Process cross sections
    if len(sections) > 0:
        with open(filename_sections, "w") as file:
            for element in sections:
                print(f"* {element}")
                _name = element.name.replace("\'", "")
                file.write(f"{_name}\n")
                file.write(f"{element.n} \t2\n")
                for i in range(element.n):
                    file.write(f"\t{element.x[i]} \t{element.y[i]} \n")
                file.write(f"\n")

    print(f"Finished writing cross-sections to '{filename_sections}'")

    return


def plot_observations(data: npt.ArrayLike, savename: str=None, keep_open: bool=False) -> plt.Figure:
    """ Visualise observation points and lines

    Input:
        data:       list of observation points and observation cross sections

    Options:
        savename:   name of figure
        keep_open:  keep figures open after finishing
    """
    ## Filename
    if not savename.endswith(".png"):
        savename += "_plot.png"

    ## Figure
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(FIGSIZE_NORMAL)
    fig.set_dpi(FIG_DPI)
    fig.set_tight_layout(True)

    ## Layout
    ax.set_title("Observation Points and Lines")
    ax.axhline(color="black", linewidth=1)
    ax.axvline(color="black", linewidth=1)

    ## Points and lines
    for element in data:
        ax = element.plot(ax)

    ## Layout part 2
    ax.grid()
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")

    ## End
    fig.savefig(savename, bbox_inches="tight", dpi=FIG_DPI, pil_kwargs=FIG_PIL_KWARGS)
    print(f"Saved figure {savename}")
    if not keep_open:
        plt.close("all")

    return fig


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath(__file__))
    obs_dir = f"{script_dir}/tests/obs"

    obs_0 = ObservationPoint(name="Center", x=0, y=0)
    obs_1 = ObservationPoint(name="Point", x=1, y=1)
    obs_2 = ObservationCrossSection(name="Line", x=[-3, 3], y=[3, 3])
    obs_3 = ObservationCrossSection(name="Diagonal", x=[-2, 3], y=[-4, 1])
    obs_4 = ObservationCrossSection(name="Curve", x=[-2, 0, -1], y=[-3, 1, 2])
    obs_5 = ObservationCrossSection(name="Square", x=[-4, -4 ,4, 4], y=[-4, 4, 4, -4])

    data = [obs_0, obs_1, obs_2, obs_3, obs_4, obs_5]

    write_observations(data, filename=obs_dir)
    plot_observations(data, savename=obs_dir)
