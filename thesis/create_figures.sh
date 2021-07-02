#!/bin/bash
# Create visualisations for all experiments

python3 ./create_visualisations.py 00 &
python3 ./create_visualisations.py 01 &
python3 ./create_visualisations.py 02 &
python3 ./create_visualisations.py 03 &
python3 ./create_visualisations.py 04 &
python3 ./create_visualisations.py 05 &
wait
echo "$(date) - Created figures exp: 00 01 02 03 04 05" >> last_runs.log

python3 ./create_visualisations.py 06 &
python3 ./create_visualisations.py 07 &
python3 ./create_visualisations.py 08 &
python3 ./create_visualisations.py 09 &
python3 ./create_visualisations.py 10 &
python3 ./create_visualisations.py 11 &
wait
echo "$(date) - Created figures exp: 06 07 08 09 10 11" >> last_runs.log

python3 ./create_visualisations.py 12 &
python3 ./create_visualisations.py 13 &
python3 ./create_visualisations.py 14 &
python3 ./create_visualisations.py 15 &
python3 ./create_visualisations.py 16 &
python3 ./create_visualisations.py 17 &
wait
echo "$(date) - Created figures exp: 12 13 14 15 16 17" >> last_runs.log

exit 0
