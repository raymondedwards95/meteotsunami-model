"""Global options and functions used in this project"""

import os
from typing import TypeAlias, Union

import dask.config
import matplotlib as mpl
import numpy as np
import seaborn as sns

# Colors
sns.set_palette(sns.color_palette("muted"))

# Matplotlib figure options
mpl.rcParams["pgf.texsystem"] = "xelatex"  # set xelatex as latex engine
mpl.rcParams["pgf.rcfonts"] = False  # unset default matplotlib fonts for pgf
mpl.rcParams["pgf.preamble"] = "\n".join(
    [
        r"\usepackage{amsmath}",
        r"\usepackage{amssymb}",
        r"\usepackage{siunitx}",
        r"\usepackage{kpfonts-otf}",
    ]
)  # apply packages for pgf

mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = "\n".join(
    [
        r"\usepackage{amsmath}",
        r"\usepackage{amssymb}",
        r"\usepackage{siunitx}",
        r"\usepackage{arev}",
    ]
)  # apply packages for png

mpl.rcParams["axes.formatter.use_mathtext"] = True  # write x10^p instead of 1ep
mpl.rcParams["axes.formatter.limits"] = -3, 4  # when to use scientific notation

# mpl.rcParams["xaxis.labellocation"] = "right"
# mpl.rcParams["yaxis.labellocation"] = "top"

mpl.rcParams["xtick.top"] = True
mpl.rcParams["xtick.bottom"] = True
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["xtick.minor.visible"] = True
mpl.rcParams["xtick.minor.size"] = 1

mpl.rcParams["ytick.left"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams["ytick.direction"] = "out"
mpl.rcParams["ytick.minor.visible"] = True
mpl.rcParams["ytick.minor.size"] = 1

mpl.rcParams["grid.alpha"] = 0.8

mpl.rcParams["legend.fancybox"] = False  # square box instead of rounded
mpl.rcParams["legend.numpoints"] = 2  # 2 points instead of 1 for plot with markers
mpl.rcParams["legend.scatterpoints"] = 1  # 1 points instead of 1 for scatter

mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["figure.titleweight"] = "bold"
mpl.rcParams["figure.titlesize"] = "x-large"

# Matplotlib video options
mpl.rcParams["animation.codec"] = "libsvtav1"
ANIM_EXT = "webm"
FFMPEG_ARGS = [
    "-auto-alt-ref",
    "0",
]

# Figure size
FIGSIZE_NORMAL = (6, 4)
FIGSIZE_LONG = (6, 9)
FIGSIZE_HIGH = (6, 6)
FIGSIZE_WIDE = (6, 3)
FIGSIZE_SMALL = (4, 3)
FIGSIZE_SQUARE = (6, 6)

# Figure save options
FIG_PIL_KWARGS = {
    "optimize": True,
    "compress_level": 9,
}

# Type Hints
Integer: TypeAlias = Union[int, np.integer]
Floating: TypeAlias = Union[float, np.floating]
Numeric: TypeAlias = Union[int, float, np.integer, np.floating]

# Paths
PATH_FUNCTIONS = os.path.dirname(os.path.realpath(__file__))
PATH_PLOTTING = f"{PATH_FUNCTIONS}/plotting"
PATH_MAIN = os.path.dirname(PATH_FUNCTIONS)
PATH_TEST = f"{PATH_MAIN}/tests"
PATH_FIGURES = f"{PATH_MAIN}/figures"
PATH_PNG = f"{PATH_MAIN}/figures/png"
PATH_PGF = f"{PATH_MAIN}/figures/pgf"

os.makedirs(PATH_TEST, exist_ok=True)
os.makedirs(PATH_FIGURES, exist_ok=True)
os.makedirs(PATH_PNG, exist_ok=True)
os.makedirs(PATH_PGF, exist_ok=True)

assert PATH_FUNCTIONS.endswith("functions")
assert PATH_PLOTTING.endswith("plotting")
assert PATH_TEST.endswith("tests")
assert PATH_FIGURES.endswith("figures")
assert PATH_PNG.endswith("png")
assert PATH_PGF.endswith("pgf")

# Dask options
dask.config.set({"array.chunk-size": "512MiB"})


# General functions
def save_figure(fig: mpl.figure, name: str, path: str = None) -> None:
    """Saves figure as a low-res png-file and a high res pgf-file

    Input:
        `fig`:      figure
        `name`:     figure name

    Options
        `path`:     path to figure
    """
    if path is None:
        path, name = os.path.split(name)

    path = path.removeprefix(PATH_MAIN)
    path = path.removeprefix("/")
    path = path.replace("/figures", "")

    print(f"## Figure path is '{path}'")
    print(f"## Figure name is '{name}'")
    print(f"## Saving figure as low resolution 'png' and full resolution 'pgf'")

    png_dir = os.path.join(PATH_PNG, path)
    pgf_dir = os.path.join(PATH_PGF, path)

    os.makedirs(png_dir, exist_ok=True)
    os.makedirs(pgf_dir, exist_ok=True)

    png_file = f"{png_dir}/{name}.png"
    pgf_file = f"{pgf_dir}/{name}.pgf"

    fig.savefig(
        png_file,
        bbox_inches="tight",
        dpi=100,
        pil_kwargs=FIG_PIL_KWARGS,
    )
    print(f"## Saved png figure as '{png_file}'")

    fig.savefig(
        pgf_file,
        bbox_inches="tight",
        dpi=1200,
    )
    print(f"## Saved pgf figure as '{pgf_file}'")

    return


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", f"{PATH_FUNCTIONS}/analysis.py"])
    subprocess.run(["python", f"{PATH_FUNCTIONS}/animation.py"])
    subprocess.run(["python", f"{PATH_FUNCTIONS}/bathymetry.py"])
    subprocess.run(["python", f"{PATH_FUNCTIONS}/observations.py"])
    subprocess.run(["python", f"{PATH_FUNCTIONS}/pressure.py"])
    subprocess.run(["python", f"{PATH_FUNCTIONS}/theory.py"])
    subprocess.run(["python", f"{PATH_FUNCTIONS}/utilities.py"])

    subprocess.run(["python", f"{PATH_PLOTTING}/base.py"])

    subprocess.run(["python", f"{PATH_PLOTTING}/alongshore.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/contour.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/crossshore.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/spectrum1d.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/spectrum2d.py"])
    subprocess.run(["python", f"{PATH_PLOTTING}/timeseries.py"])
