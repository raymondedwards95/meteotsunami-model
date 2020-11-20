# Reproduction of experiments in paper An et al. (2012)

Scripts are used to create certain fields:
* `create_bathymetry.py` creates the `bathymetry.xyb` file;
* `create_pressure.py` sets up a space and time varying pressure field in `pressure.amp`;
* `create_grid.py` creates a grid in `FlowFM_net.nc` using macros in `Delta Shell` (note that `Delta Shell` should be used to run this file).

Other files are `forcing.ext` which activates the pressure field in `pressure.amp` and `FlowFM.mdu` is the main file with parameters for simulations. Normally these two files should not change much.
