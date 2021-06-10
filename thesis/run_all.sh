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

echo "Last runs:" > last_runs.log

# run cases
bash ./run_exp.sh 00

exit 0
