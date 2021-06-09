""" Creates a file that contains the observation cross-sections for D3D-FM-FLOW """

import os
import sys

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.observations as fo


## Parameters
x_min = 0
x_max = 1e6

y_min = -1e7
y_max = +1e7

obs_dir = f"{os.path.dirname(os.path.realpath(__file__))}/obs"
os.makedirs(obs_dir, exist_ok=True)


## Observation Points and Cross Sections
obs = []
obs.append(
    fo.ObservationPoint("Test Point", x=10e3, y=x_max/2.)
)

for i in range(4):
    _x = i * 1e4
    obs.append(
        fo.ObservationCrossSection(
            name=f"Along Shore x={_x/1000 :0.0f}km",
            x=[_x, _x],
            y=[y_min, y_max]
        )
    )

for i in range(5):
    _y = i*y_max/4
    obs.append(
        fo.ObservationCrossSection(
            name=f"Cross Shore y={_y/1000:0.0f}km",
            x=[x_min, x_max],
            y=[_y, _y]
        )
    )

fo.write_observations(obs, filename=f"{obs_dir}/exp_00")
print("Finished creating observation-files")
