#!/bin/bash
# Runs and make figures for all cases with name 'input_exp_**.mdu' in this folder
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
bash ./run_exp.sh 00 &
bash ./run_exp.sh 01 &
bash ./run_exp.sh 02 &
bash ./run_exp.sh 03 &
bash ./run_exp.sh 04 &
bash ./run_exp.sh 05 &
wait

bash ./run_exp.sh 06 &
bash ./run_exp.sh 07 &
bash ./run_exp.sh 08 &
bash ./run_exp.sh 09 &
bash ./run_exp.sh 10 &
bash ./run_exp.sh 11 &
wait

bash ./run_exp.sh 12 &
bash ./run_exp.sh 13 &
bash ./run_exp.sh 14 &
bash ./run_exp.sh 15 &
bash ./run_exp.sh 16 &
bash ./run_exp.sh 17 &
wait

echo "$(date) - Finished all runs!" >> last_runs.log

exit 0
