# Reproduction of experiments in paper An et al. (2012)

This part is aimed at reproducing the results of the experiments in the paper by An et al. (2012).
The purpose is to get the model running and to get results comparable with those in the mentioned paper.
An other goal is to investigate the effects of the different model parameters on the simulations.


### Structure

Scripts are used to create all different fields:
* `create_bathymetry.py` creates `*.xyb` (and additional `*.nc`) files in `./bathymetry/` that contain the bathymetry;
* `create_grid.py` creates all grids in `./grid/` using macros in `Delta Shell` (this also means that `Delta Shell` is necessary to run this file);
* `create_observations.py` creates `*.xyn` and `*_crs.pli` in `./obs/` that contain points and cross-sections for data aquisition during simulation;
* `create_pressure.py` sets up a space and time varying pressure field and writes them to `.amp` (and additionally in `*.nc`) files in `./pressure/` and also creates corresponding `forcing_*.ext` files in that location to activate the pressure fields in the model.

The `input*.mdu` are the main files that contain all parameters for simulations. They also specify what additional files are used, like bathymetry, pressure and grid.

Three scripts are used to visualize data:
* `create_animations.py` for animations;
* `create_comparison.py` for comparisons with the paper by An et al. (2012);
* `create_figures.py` for figures.


### Running experiments

Simulations can be done by running `bash run_reproduction.sh case`, where `case` is a case number as in the table below.
It will run the model and afterwards it will create some figures.

Running `bash run_all.sh` will run all experiments.


### Experiments

Multiple simulations are done, related to the base experiments.
In these cases one or more parameters are changed, compared to the base experiment.
All different cases are numbered.

| Number | Changes | Comments |
| :--- | :--- | :--- |
||||
| 00 | - | Base experiment |
||||
| 01 | Change `Dtmax` from `20` to `10` s | Maximum computational time step |
| 02 | Change `Dtmax` from `20` to `60` s | See 01 |
||||
| 03 | Change `Teta0` from `0.55` to `0.99` | Parameter for implicitness of the numerical scheme (between 0.5 and 1.0) |
| 04 | Change `Teta0` from `0.55` to `0.59` | See 03 |
| 05 | Change `Teta0` from `0.55` to `0.51` | See 03 |
||||
| 10 | Set `dx_p` and `dy_p` from `10` to `5` km | Spatial resolution of pressure distribution |
| 11 | Set `dx_p` and `dy_p` from `10` to `20` km | See 10 |
| 12 | Set `dx_p` and `dy_p` from `10` to `40` km | See 10 |
||||
| 15 | Set `dt_p` from `30` to `60` min | Temporal resolution of pressure distribution |
| 16 | Set `dt_p` from `30` to `20` min | See 15 |
| 17 | Set `dt_p` from `30` to `10` min | See 15 |
||||
| 20 | Set `dx` and `dy` from `10` to `20` km | Spatial resolution of model |
| 21 | Set `dx` and `dy` from `10` to `40` km | See 20 |
| 22 | Set `dx` and `dy` from `10` to `5` km | See 20, skipped |
||||
| 31 | Set `x0_p` from `0` km to `100` km | Position of center of pressure disturbance |
| 32 | Set `x0_p` from `0` km to `200` km | See 31 |
| 33 | Set `x0_p` from `0` km to `300` km | See 31 |
||||
| 36 | Set `alpha` from `1/400` to `1/800`  | Slope of bottom |
| 37 | Set `alpha` from `1/400` to `1/800` and set `x0_p` from `0` km to `100` km | See 31 and 36 |
| 38 | Set `alpha` from `1/400` to `1/800` and set `x0_p` from `0` km to `200` km | See 31 and 36 |
| 39 | Set `alpha` from `1/400` to `1/800` and set `x0_p` from `0` km to `300` km | See 31 and 36 |
||||
| 41 | Set `alpha` from `1/400` to `0` and set `average_depth` to `250` m and set `x0_p` from `0` to `500` km | Flat bottom and see 33|
| 42 | Set `alpha` from `1/400` to `0` and set `average_depth` to `100` m and set `x0_p` from `0` to `500` km | See 41 |
| 43 | Set `alpha` from `1/400` to `0` and set `average_depth` to `500` m and set `x0_p` from `0` to `500` km | See 41 |
