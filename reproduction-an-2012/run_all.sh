#!/bin/bash
# Runs and make figures for all cases with name 'input_repr_**.mdu' in this folder
# make sure the settings in the root folder are correct, following the template

# set path to model executables
if [ ! -f ../settings.sh ]
  then
    echo "../settings.sh does not exist. Creating file... Probably it needs to be adjusted manually."
    cp ../setup/settings.sh.template ../settings.sh
fi

source ../settings.sh
echo "Path to executables is $DFLOWFM_BIN_PATH"

# check if executable exists
if [ ! -f $INPUT_FILE ]
  then
    echo "$DFLOWFM_BIN_PATH/run_dflowfm.sh does not exist. Check ../setup/settings.sh for errors. Exiting..."
    exit 1
fi

echo "$(date) - Starting new runs!" >> last_runs.log

# run cases
bash ./run_reproduction.sh 00 &
bash ./run_reproduction.sh 01 &
bash ./run_reproduction.sh 02 &
bash ./run_reproduction.sh 03 &
bash ./run_reproduction.sh 04 &
bash ./run_reproduction.sh 05 &
wait

bash ./run_reproduction.sh 10 &
bash ./run_reproduction.sh 11 &
bash ./run_reproduction.sh 12 &
bash ./run_reproduction.sh 15 &
bash ./run_reproduction.sh 16 &
bash ./run_reproduction.sh 17 &
wait

bash ./run_reproduction.sh 31 &
bash ./run_reproduction.sh 32 &
bash ./run_reproduction.sh 33 &
bash ./run_reproduction.sh 36 &
bash ./run_reproduction.sh 37 &
bash ./run_reproduction.sh 38 &
bash ./run_reproduction.sh 39 &
wait

bash ./run_reproduction.sh 41 &
bash ./run_reproduction.sh 42 &
bash ./run_reproduction.sh 43 &

bash ./run_reproduction.sh 20 &
bash ./run_reproduction.sh 21 &
# bash ./run_reproduction.sh 22 &
wait

python3 ./create_comparison.py

exit 0
