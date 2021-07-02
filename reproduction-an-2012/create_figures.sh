#!/bin/bash
# Create visualisations for all reproduction-experiments

echo "$(date) - Starting creation of figures repr"
echo "$(date) - Starting creation of figures repr" >> last_runs.log

python3 ./create_visualisations.py 00 &
python3 ./create_visualisations.py 01 &
python3 ./create_visualisations.py 02 &
python3 ./create_visualisations.py 03 &
python3 ./create_visualisations.py 04 &
python3 ./create_visualisations.py 05 &
wait
echo "$(date) - Created figures repr: 00 01 02 03 04 05" 
echo "$(date) - Created figures repr: 00 01 02 03 04 05" >> last_runs.log

python3 ./create_visualisations.py 10 &
python3 ./create_visualisations.py 11 &
python3 ./create_visualisations.py 12 &
python3 ./create_visualisations.py 15 &
python3 ./create_visualisations.py 16 &
python3 ./create_visualisations.py 17 &
wait
echo "$(date) - Created figures repr: 10 11 12 15 16 17" 
echo "$(date) - Created figures repr: 10 11 12 15 16 17" >> last_runs.log

python3 ./create_visualisations.py 20 &
python3 ./create_visualisations.py 21 &
python3 ./create_visualisations.py 31 &
python3 ./create_visualisations.py 32 &
python3 ./create_visualisations.py 33 &
python3 ./create_visualisations.py 36 &
wait
echo "$(date) - Created figures repr: 20 21 31 32 33 36" 
echo "$(date) - Created figures repr: 20 21 31 32 33 36" >> last_runs.log

python3 ./create_visualisations.py 37 &
python3 ./create_visualisations.py 38 &
python3 ./create_visualisations.py 39 &
python3 ./create_visualisations.py 41 &
python3 ./create_visualisations.py 42 &
python3 ./create_visualisations.py 43 &
wait
echo "$(date) - Created figures repr: 37 38 39 41 42 43" 
echo "$(date) - Created figures repr: 37 38 39 41 42 43" >> last_runs.log

wait
echo "$(date) - Created all figures repr"
echo "$(date) - Created all figures repr" >> last_runs.log

exit 0
