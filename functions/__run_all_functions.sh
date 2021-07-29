#!/bin/bash
# Run the 'main' parts of the files containing functions for testing purposes

python3 ./__init__.py
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
