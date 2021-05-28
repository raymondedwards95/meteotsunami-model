""" Script to compare different simulations """

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import os
import sys

# fix for importing functions below
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import functions.visualisation as fv


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
assert len(titles) == num_comparisons


### Paths
file_dir = os.path.dirname(os.path.realpath(__file__))
output_dir = file_dir + "/output/"
figure_dir = file_dir + "/figures/"


### Defining Visualisations
def vis_alongshore(data_list, title):
    raise NotImplementedError


def vis_crossshore(data_list, title):
    raise NotImplementedError


def make_comparison(cases, title):
    raise NotImplementedError


def make_comparison(cases, title, id="test"):
    # Paths
    if isinstance(id, float):
        id = f"{id:03.0f}"
    savename = figure_dir + f"/comparison_{id}"
    os.makedirs(savename, exist_ok=True)

    # Import data (lazy?)
    data_list = []
    for i in range(len(cases)):
        case = f"{cases[i]:02.0f}"
        print(f"Reading data for case {case}")
        data_list.append(xr.open_dataset(f"{output_dir}/data_repr_{case}.nc"))
    
    # Make figures
    vis_alongshore(data_list, title, cases, savename)
    vis_crossshore(data_list, title, cases, savename)

    return


### Executing Comparison
for i in range(num_comparisons):
    make_comparison(cases[i], titles[i])
