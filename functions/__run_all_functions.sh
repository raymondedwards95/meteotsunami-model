#!/bin/bash
# Run the 'main' parts of the files containing functions for testing purposes

echo "Running all python-files in this directory $PWD"

python3 ./analysis.py
python3 ./animation.py
python3 ./bathymetry.py
python3 ./observations.py
python3 ./pressure.py
python3 ./regrid.py
python3 ./utilities.py
python3 ./visualisation.py
wait

exit 0
