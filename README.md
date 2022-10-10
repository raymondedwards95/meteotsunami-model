# meteotsunami-model

Scripts and files to run simulations in Delft3d-FM.


### Structure

The main parts are the folders `reproduction-an-2012` and `thesis` as they contain the input files, output files and results for all different simulations.
The `reproduction-an-2012` is aimed at setting up model parameters.
The main research part is in `thesis`, which is more focused at physical parameters.

The folder `functions` contain some methods to prepare, process and visualise simulations, while `setup` contains some other useful templates.


### Remarks

To run an experiment with D3D-FM, use for example:
```
/path/to/d3d-fm/lnx64/bin/run_dflowfm.sh /path/to/mdufile.mdu
```
