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
mpl.rcParams["animation.codec"] = "libaom-av1"
ANIM_EXT = "webm"

# Figure size
FIGSIZE_NORMAL = (7, 4)
FIGSIZE_LONG = (7, 9)
FIGSIZE_HIGH = (7, 6)
FIGSIZE_WIDE = (7, 3)
FIGSIZE_SMALL = (3, 2)
FIGSIZE_SQUARE = (7, 7)

FIG_DPI = 300

# Figure save options
FIG_PIL_KWARGS = {"optimize": True, "compress_level": 9}

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
