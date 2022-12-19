""" Creates a file that contains the observation cross-sections for D3D-FM-FLOW """

import os
import sys

# fmt: off
# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from functions import *
import functions.observations as fo
# fmt: on


# Parameters
x_min = 0.0
x_max = 1e7

y_min = -1e7
y_max = +1e7

script_dir = os.path.dirname(os.path.realpath(__file__))
obs_dir = f"{script_dir}/obs"
os.makedirs(obs_dir, exist_ok=True)


# Observation Points and Cross Sections
obs = []
obs.append(
    fo.ObservationPoint(
        name="Test Point",
        x=10e3,
        y=x_max / 2.0,
    )
)

for i in range(4):
    _x = i * 1e4
    obs.append(
        fo.ObservationCrossSection(
            name=f"Along: $x={_x/1e3 :0.0f}$ km",
            x=[_x, _x],
            y=[y_min, y_max],
        )
    )

for i in range(5):
    _y = i * y_max / 4
    obs.append(
        fo.ObservationCrossSection(
            name=f"Cross: $y={_y/1e6:0.1f}$ Mm",
            x=[x_min, x_max],
            y=[_y, _y],
        )
    )

fo.write_observations(
    obs,
    filename=f"{obs_dir}/repr_00",
)
fo.plot_observations(
    obs,
    savename=f"{obs_dir}/repr_00",
    scale="Mm",
)
print("Finished creating observation-files")
