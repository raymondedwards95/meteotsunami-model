#!/bin/bash
# Run all scripts for creating necessary files

python3 ./create_bathymetry.py &
python3 ./create_observations.py &
python3 ./create_pressure.py &

exit 0
