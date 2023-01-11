# meteotsunami-model

Scripts and files to run simulations in Delft3d-FM.


### Structure

The main parts are the folders `reproduction-an-2012` and `thesis` as they contain the input files, output files and results for all different simulations.
The `reproduction-an-2012` is aimed at setting up model parameters.
The main research part is in `thesis`, which is more focused at physical parameters.

In the folder `theory` some scripts are included to analyse theoretical relations between different parameters and variables.

The folder `functions` contain methods to prepare, process and visualise simulations, while `setup` contains some other useful templates.


### Setup

Running `main_script.sh` once will copy setup files into the root of this folder. This new `settings.sh` should include the path to the executables.


### Running

To run all experiments and to create all visualisations, `main_script.sh --all` can be used. Running `main_script.sh --help` shows information about the use of the script.


### Other remarks

To run a single experiment with D3D-FM, use for example:
```
/path/to/d3d-fm/lnx64/bin/run_dflowfm.sh /path/to/mdufile.mdu
```


### Versions used

The following versions of software are used. The project may run on other versions, but that might give different results or even errors.

* Delft3D-FM 0.9.1 (8 July 2020);
* Python 3.11;

and the packages in `requirements.txt`:
