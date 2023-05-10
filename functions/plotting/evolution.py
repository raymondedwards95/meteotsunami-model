import os
import sys

import matplotlib.pyplot as plt
import xarray as xr

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
from functions import *
# fmt: on


script_dir = f"{PATH_PLOTTING}"
figure_dir = f"{PATH_TEST}"

data_a = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_00.nc"
data_b = f"{PATH_MAIN}/reproduction-an-2012/output/data_repr_17.nc"

data_a = xr.open_dataset(data_a, chunks="auto")
data_b = xr.open_dataset(data_b, chunks="auto")

x_ref = 1e4

fig, axes = plt.subplots(2, 1, sharex=True)

time_max = 0

for i, data in enumerate([data_a, data_b]):
    for j, var in enumerate(["wl", "p"]):
        data_time = data["t"].values.astype("datetime64[s]").astype(float) / 3600.0
        data_max = data[var].interp(x=x_ref).max(dim=["y"])
        data_min = -1.0 * data[var].interp(x=x_ref).min(dim=["y"])

        a = data_max > 0.99 * data_max.max()
        b = data_min > 0.99 * data_min.max()

        time_max = np.max(
            [
                time_max,
                data_time[a.argmax().values],
                data_time[b.argmax().values],
            ]
        )

        axes[j].plot(
            data_time,
            data_max,
            label=f"{i} max",
        )
        axes[j].fill_between(
            data_time,
            data_max,
            alpha=0.02,
        )

        axes[j].plot(
            data_time,
            data_min,
            label=f"{i} min",
        )
        axes[j].fill_between(
            data_time,
            data_min,
            alpha=0.02,
        )

for ax in axes:
    ax.set_xlim(0, time_max)
    ax.set_ylim(0, None)
    ax.legend()
    ax.grid()

axes[1].set_xlabel("Time [hours]")
axes[0].set_ylabel("Water Level \n[meters]")
axes[1].set_ylabel("Surface Air Pressure \n[pascals]")

fig.align_labels()

plt.savefig(
    f"{figure_dir}/test_evolution",
    bbox_inches="tight",
    pil_kwargs=FIG_PIL_KWARGS,
)
