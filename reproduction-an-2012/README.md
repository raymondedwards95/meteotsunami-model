# Reproduction of experiments in paper An et al. (2012)

Scripts are used to create certain fields:
* `create_bathymetry.py` creates the `bath*.xyb` file;
* `create_pressure.py` sets up a space and time varying pressure field in `pressure*.amp`;
* `create_grid.py` creates a grid in `grid*.nc` using macros in `Delta Shell` (note that `Delta Shell` should be used to run this file). Furthermore it also creates observation points and observation cross sections in the files `obs*.xyn` and `obs*_crs.pli`.

Other files are `forcing*.ext` which activates the pressure field in `pressure*.amp`, and `input*.mdu` is the main file with parameters for simulations. Normally these two files should not change much.

### Experiments

Multiple simulations are done, related to the base experiments. In these cases one or more parameters are changed, compared to the base experiment. All different cases are numbered.

| Number | Changes | Comments |
| :--- | :--- | :--- |
| 00 |  | Base experiment |
| 01 | Change `Dtmax` from `20` to `10` s | Smaller maximum compuational time step |
| 02 | Change `Dtmax` from `20` to `3` s | See 01 |
| 03 | Change `Teta0` from `0.55` to `0.99` | Parameter for implicitness of the numerical scheme (between 0.5 and 1.0) |
| 04 | Change `Teta0` from `0.55` to `0.59` | See 03 |
| 05 | Change `Teta0` from `0.55` to `0.51` | See 03 |
