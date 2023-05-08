import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
# fmt: on


script_dir = f"{PATH_PLOTTING}"
figure_dir = f"{PATH_TEST}"

data_a = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_00.nc"
data_b = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_01.nc"

data_a = xr.open_dataset(data_a, chunks="auto")
data_b = xr.open_dataset(data_b, chunks="auto")


plt.figure()
plt.plot(
    data_a["t"].values.astype("datetime64[s]").astype(float) / 3600.0,
    data_a["wl"].max(dim=["x", "y"]),
    label="00 max",
)
plt.plot(
    data_a["t"].values.astype("datetime64[s]").astype(float) / 3600.0,
    -1.0 * data_a["wl"].min(dim=["x", "y"]),
    label="00 min",
)
plt.plot(
    data_b["t"].values.astype("datetime64[s]").astype(float) / 3600.0,
    data_b["wl"].max(dim=["x", "y"]),
    label="01 max",
)
plt.plot(
    data_b["t"].values.astype("datetime64[s]").astype(float) / 3600.0,
    -1.0 * data_b["wl"].min(dim=["x", "y"]),
    label="01 min",
)
plt.ylim(0, None)
plt.legend()
plt.xlabel("Time [hours]")
plt.ylabel("Water Level [meters]")
plt.savefig(f"{figure_dir}/test_evolution")
