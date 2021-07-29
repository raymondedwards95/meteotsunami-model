#!/bin/bash
# Run the 'main' parts of the files containing functions for testing purposes

echo "Running all python-files in this directory $PWD"

echo ">> analysis.py"
python3 ./analysis.py

echo ">> animation.py"
python3 ./animation.py

echo ">> bathymetry.py"
python3 ./bathymetry.py

echo ">> observations.py"
python3 ./observations.py

echo ">> pressure.py"
python3 ./pressure.py

echo ">> regrid.py"
python3 ./regrid.py

echo ">> utilities.py"
python3 ./utilities.py

echo ">> visualisation.py"
python3 ./visualisation.py

wait
exit 0
