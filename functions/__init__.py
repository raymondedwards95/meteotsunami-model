""" Global options """
import os
from typing import TypeAlias, Union

import dask.config
import matplotlib as mpl
import numpy as np
import seaborn as sns

# Colors
sns.set_palette(sns.color_palette("muted"))

# Matplotlib figure options


# Matplotlib video options
mpl.rcParams["animation.codec"] = "vp8"
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

FIG_DPI = 300

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
        dpi=75,
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
