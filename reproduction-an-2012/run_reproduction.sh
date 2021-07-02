#!/bin/bash
# Runs and make figures for a given case, for which the path to the inputfile (*.mdu) should be given as an argument
# Make sure the settings in the root folder are correct, following the template
# Also check if required files are present; use the python-scripts if necessary

# check if there are arguments
if [ $# -eq 0 ]
  then
    echo "No input file given. Exiting..."
    exit 1
fi

# set case number and input file name
CASE_NUMBER=$1
INPUT_FILE="./input_repr_$CASE_NUMBER.mdu"
echo "Input case is $CASE_NUMBER"
echo "Input file is $INPUT_FILE"

# check if input file exists
if [ ! -f $INPUT_FILE ]
  then
    echo "Input file $INPUT_FILE does not exist. Exiting..."
    exit 1
fi

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

# start timer
T_START=$(date +%s)
echo "$(date) - Start simulation repr case $CASE_NUMBER "
echo "$(date) - Start simulation repr case $CASE_NUMBER " >> last_runs.log

# run simulation
$DFLOWFM_BIN_PATH/run_dflowfm.sh $INPUT_FILE

# process outputs
python3 create_visualisations.py $CASE_NUMBER --reprocess 1

# stop timer
T_END=$(date +%s)
let TIME_DELTA="$T_END-$T_START"
let TIME_DELTA_MINUTE="$TIME_DELTA/60"
let TIME_DELTA_HOUR="$TIME_DELTA/3600"

# end script
echo "Finished simulations and processing and visualising data in $TIME_DELTA seconds ($TIME_DELTA_MINUTE minutes or $TIME_DELTA_HOUR hours)"
echo "$(date) - End simulation repr case $CASE_NUMBER "
echo "$(date) - End simulation repr case $CASE_NUMBER " >> last_runs.log
exit 0
