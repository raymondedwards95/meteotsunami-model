""" Creates a file that contains a grid for D3D-FM-FLOW 
Note: this script should be run from Delta-Shell
"""

import os
import sys

import Libraries.FlowFlexibleMeshFunctions as FFMF
import Libraries.StandardFunctions as SF


## Grid
x_min = 0
x_max = 1e6
x_step = 1e4

y_min = -1e7
y_max = +1e7
y_step = 1e4


## Write to file
x_length = x_max - x_min
x_num = int(x_length / x_step) + 1

y_length = y_max - y_min
y_num = int(y_length / y_step) + 1

model = FFMF.WaterFlowFMModel()
SF.AddToProject(model)
FFMF.GenerateRegularGridForModelUsingExtend(
    model, x_length, y_length, x_num, y_num, x_min, y_min)
model.WriteNetFile(os.path.dirname(os.path.realpath(__file__)) + "/FlowFM_net.nc")


## Observation Points and Cross Sections
obs = []
obs.append(fo.ObservationPoint("Test Point", x=10e3, y=x_max/2.))

for i in range(4):
    _x = i * 1e4
    obs.append(fo.ObservationCrossSection(name=f"Along Shore x={_x/1000:0.0f}km", x=[_x, _x], y=[y_min, y_max]))

for i in range(5):
    _y = i*y_max/4
    obs.append(fo.ObservationCrossSection(name=f"Cross Shore y={_y/1000:0.0f}km", x=[x_min, x_max], y=[_y, _y]))

## Write to file
for i in range(num_cases):
    x_length = x_max - x_min
    x_num = int(x_length / x_step[i]) + 1

    y_length = y_max - y_min
    y_num = int(y_length / y_step[i]) + 1

    model = FFMF.WaterFlowFMModel()
    SF.AddToProject(model)
    FFMF.GenerateRegularGridForModelUsingExtend(
        model, x_length, y_length, x_num, y_num, x_min, y_min)
    model.WriteNetFile(os.path.dirname(os.path.realpath(__file__)) + "/grid_repr_{:2.0f}_net.nc".format(case[i]))
