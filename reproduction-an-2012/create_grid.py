""" Creates a file that contains a grid for D3D-FM-FLOW
Note: this script should be run from Delta-Shell
Note: this script does not work in the console version of Delta-Shell
best solution for now is to copy the contents of this file into the Delta-Shell gui
"""

import os

import Libraries.FlowFlexibleMeshFunctions as FFMF
import Libraries.StandardFunctions as SF


## Filename (adjust following!)
grid_dir = os.path.dirname(os.path.realpath(__file__)) + "/grid"
os.makedirs(grid_dir, exist_ok=True)


## Grid
case = [0, 20, 21, 22]
num_cases = len(case)

x_min = 0
x_max = 1e6
x_step = [1e4, 2e4, 4e4, 5000]

y_min = -1e7
y_max = +1e7
y_step = x_step

assert len(x_step) == num_cases
assert len(y_step) == num_cases

x_length = x_max - x_min
y_length = y_max - y_min


## Write to file
for i in range(num_cases):
    x_num = int(x_length / x_step[i]) + 1
    y_num = int(y_length / y_step[i]) + 1

    net_file_path = grid_dir + "/repr_{:02.0f}_net.nc".format(float(case[i]))

    model = FFMF.WaterFlowFMModel()
    SF.AddToProject(model)
    model.WriteNetFile(net_file_path)
    FFMF.GenerateRegularGridForModelUsingExtend(model, x_length, y_length, x_num, y_num, x_min, y_min)
    model.WriteNetFile(net_file_path)

    print("Finished writing to", net_file_path)

print("Finished writing grid")
