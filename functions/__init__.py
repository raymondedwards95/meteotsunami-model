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
FIGSIZE_NORMAL = (7, 4)
FIGSIZE_LONG = (7, 9)
FIGSIZE_HIGH = (7, 6)
FIGSIZE_WIDE = (7, 3)
FIGSIZE_SMALL = (3, 2)
FIGSIZE_SQUARE = (7, 7)

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
functions_dir = os.path.dirname(os.path.realpath(__file__))
main_dir = os.path.dirname(functions_dir)
test_results_dir = f"{functions_dir}/tests"
os.makedirs(test_results_dir, exist_ok=True)

assert functions_dir.endswith("functions")
assert test_results_dir.endswith("tests")

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

    print(f"## Figure path is '{path}'")
    print(f"## Figure name is '{name}'")
    print(f"## Saving figure as low resolution 'png' and full resolution 'pgf'")

    png_dir = os.path.join(path, "png")
    pgf_dir = os.path.join(path, "pgf")

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
