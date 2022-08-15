""" Global options """
from typing import TypeAlias, Union
import numpy as np
import seaborn as sns

### Colors
sns.set_palette(sns.color_palette("muted"))

### Figure size
FIGSIZE_NORMAL = (7, 4)
FIGSIZE_LONG = (7, 9)
FIGSIZE_HIGH = (7, 6)
FIGSIZE_WIDE = (7, 3)
FIGSIZE_SMALL = (3, 2)
FIGSIZE_SQUARE = (7, 7)

FIG_DPI = 300

### Type Hints
Numeric: TypeAlias = Union[int, float, np.integer, np.floating]
