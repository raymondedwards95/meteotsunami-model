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
