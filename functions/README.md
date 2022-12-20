# Functions

This folder contains functions to help setting up simulations, to process model outputs and to visualise model output.

* `analysis.py` provides some tools to analyse data and to compute certain parameters;
* `animation.py` provides methods to create animations of model output;
* `bathymetry.py` helps with setting up the bathymetry;
* `observations.py` contains functions to write files with observation points and cross sections;
* `pressure.py` helps with setting up space and time varying pressure fields;
* `regrid.py` contains functions to convert D3D-FM output from a flexible mesh to a fixed regular grid (it can also be used as standalone);
* `theory.py` contains theoretical relationships;
* `utilities.py` provides some useful general functions;

Importing this folder, i.e. `__init__.py` provides some general parameters and options.

Functions for visualising are contained in the folder `plotting`. These are

* `plotting/contour.py` for visualising spatial variations;
* `plotting/profiles.py` for graphing profiles;
* `plotting/spectrum.py` for computing and visulalising 1d and 2d power spectra;
* `plotting/timeseries.py` for visualising temporal variations.
