""" Script to compare different simulations """

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


### Cases
cases = [
    [0, 1, 2],
    [0, 3, 4, 5],
    [0, 10, 11, 12],
    [0, 15, 16, 17],
    [0, 20, 21],
]

titles = [
    "Maximum Computational Time Step",
    "Implicitness of Numerical Scheme",
    "Spatial Resolution - Pressure",
    "Temporal Resolution - Pressure",
    "Spatial Resolution - Model"
]

num_comparisons = len(cases)


### Defining Visualisations
def vis_alongshore(data_list, title):
    raise NotImplementedError


def vis_crossshore(data_list, title):
    raise NotImplementedError


def make_comparison(cases, title):
    raise NotImplementedError


### Executing Comparison
for i in range(num_comparisons):
    make_comparison(cases[i], titles[i])
