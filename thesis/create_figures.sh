
# run cases
python3 ./create_visualisation.py 00 &
python3 ./create_visualisation.py 01 &
python3 ./create_visualisation.py 02 &
python3 ./create_visualisation.py 03 &
python3 ./create_visualisation.py 04 &
python3 ./create_visualisation.py 05 &
wait
echo "Created figures exp: 00 01 02 03 04 05\n" >> last_runs.log

python3 ./create_visualisation.py 06 &
python3 ./create_visualisation.py 07 &
python3 ./create_visualisation.py 08 &
python3 ./create_visualisation.py 09 &
python3 ./create_visualisation.py 10 &
python3 ./create_visualisation.py 11 &
wait
echo "Created figures exp: 06 07 08 09 10 11\n" >> last_runs.log

python3 ./create_visualisation.py 12 &
python3 ./create_visualisation.py 13 &
python3 ./create_visualisation.py 14 &
python3 ./create_visualisation.py 15 &
python3 ./create_visualisation.py 16 &
python3 ./create_visualisation.py 17 &
wait
echo "Created figures exp: 12 13 14 15 16 17\n" >> last_runs.log

exit 0
