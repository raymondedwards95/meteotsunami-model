import os

import Libraries.FlowFlexibleMeshFunctions as FFMF
import Libraries.StandardFunctions as SF


x_min = 0
x_max = 1e6
x_step = 1e4

y_min = -1e7
y_max = +1e7
y_step = 1e4


x_length = x_max - x_min
x_num = int(x_length / x_step) + 1

y_length = y_max - y_min
y_num = int(y_length / y_step) + 1

model = FFMF.WaterFlowFMModel()
SF.AddToProject(model)
FFMF.GenerateRegularGridForModelUsingExtend(
    model, x_length, y_length, x_num, y_num, x_min, y_min)
model.WriteNetFile(os.path.dirname(os.path.realpath(__file__)) + "/FlowFM_net.nc")
