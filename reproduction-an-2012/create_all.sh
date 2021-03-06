#!/bin/bash
# Run all scripts for creating necessary files

echo "$(date) - Creating additional input files for repr"
echo "$(date) - Creating additional input files for repr" >> last_runs.log

# Run scripts
python3 ./create_bathymetry.py
python3 ./create_observations.py
python3 ./create_pressure.py
wait
echo "$(date) - Finished creating additional input files for repr"
echo "$(date) - Finished creating additional input files for repr" >> last_runs.log

exit 0
