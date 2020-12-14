#!/bin/bash
# Runs and make figures for a given case, for which the path to the inputfile (*.mdu) should be given as an argument
# 
# make sure the settings in the root folder are correct, following the template

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
source ../setup/settings.sh
echo "Path to executables is $DFLOWFM_BIN_PATH"

# check if executable exists
if [ ! -f $INPUT_FILE ]
  then
    echo "$DFLOWFM_BIN_PATH/run_dflowfm.sh does not exist. Check ../setup/settings.sh for errors. Exiting..."
    exit 1
fi

# start timer
T_START=$(date +%s)

# run simulation
$DFLOWFM_BIN_PATH/run_dflowfm.sh $INPUT_FILE

# process outputs
python3 create_visualisations.py $CASE_NUMBER

# stop timer
T_END=$(date +%s)
let TIME_DELTA="$T_END-$T_START"
let TIME_DELTA_HOUR="$TIME_DELTA/3600"

# end script
echo "Finished simulations and processing and visualising data in $TIME_DELTA seconds ($TIME_DELTA_HOUR hours)"
exit 0
