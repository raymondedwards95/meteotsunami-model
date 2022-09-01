#!/bin/bash
# Runs all simulations in the given folders.
# Also processes the model output.

# Maximum number of simultanious tasks
Ntasks=5

# List of folders with inputfiles
FolderList=("./reproduction-an-2012 ./thesis")

# Global settings
if [ ! -f ./settings.sh ]
  then
    echo "./settings.sh does not exist. Creating file... Probably it needs to be adjusted manually."
    cp ./setup/settings.sh.template ./settings.sh
    exit 1
fi

BaseDir="${PWD}"
LogFile="${PWD}/last_runs.log"
LogFolder="${PWD}/logs"

mkdir -p "${LogFolder}"

echo
echo "### Prepare simulations"
echo "## Reading './settings.sh'"
source ./settings.sh
echo "# Current directory is '$BaseDir'"
echo "# Logging file is '$LogFile'"
echo "# Path to executables is $DFLOWFM_BIN_PATH"
echo "# Maximum number of simultanious tasks is $Ntasks"

if [ ! -f "$DFLOWFM_BIN_PATH/run_dflowfm.sh" ]
  then
    echo "$DFLOWFM_BIN_PATH/run_dflowfm.sh does not exist. Check ./setup/settings.sh. Exiting..."
    exit 1
fi

# Define computations
function func_computations ()
{
    # Catch arguments
    local LocalCase=$1
    local LocalInputFile=$2
    local LocalOutputFile=$3
    local LocalRegridFile=$4
    local LocalIdentifier=$5

    # Start timer
    local TStart=$(date +%s)

    echo "  '$LocalInputFile' ->"
    echo "  '$LocalOutputFile' ->"
    echo "  '$LocalRegridFile'"

    # Do computations
    echo "# Start simulations for $LocalIdentifier case $LocalCase"
    echo "$(date) - Start simulation for '$LocalInputFile'" >> $LogFile
    $DFLOWFM_BIN_PATH/run_dflowfm.sh $LocalInputFile 1> "${LogFolder}/sim_${LocalIdentifier}_${LocalCase}.log" 2>&1

    local TMiddle=$(date +%s)

    # Regrid data
    echo "# Start Regridding for $LocalIdentifier case $LocalCase"
    echo "$(date) - Start regridding for '$LocalInputFile'" >> $LogFile
    python3 ${BaseDir}/functions/regrid.py $LocalOutputFile $LocalRegridFile --delete-original-model-output 1 1> "${LogFolder}/regrid_${LocalIdentifier}_${LocalCase}.log" 2>&1

    # End
    local TEnd=$(date +%s)
    local dTComputation=$(python3 -c "print(f'{($TMiddle - $TStart) / 60:0.1f}')")
    local dTRegrid=$(python3 -c "print(f'{($TEnd - $TMiddle) / 60:0.1f}')")
    local dTTotal=$(python3 -c "print(f'{($TEnd - $TStart) / 60:0.1f}')")

    echo "# Finished computations for $LocalIdentifier case $LocalCase"
    echo "# Case $LocalIdentifier $LocalCase took $dTTotal minutes in total: $dTComputation minutes for D3D and $dTRegrid minutes for regridding"
    echo "$(date) - Finished all for '$LocalInputFile'" >> $LogFile
}

# Loop over all folders
for Folder in $FolderList
do
    # Move to folder
    cd $Folder
    echo
    echo "### Doing simulations in $PWD"

    # Find cases
    FileName=$(find *.mdu | grep -m 1 -P -o "^\D+")
    Identifier=$(find *.mdu | grep -m 1 -P -o "_(\w+)_" | grep -P -o "[^_]*")
    CaseNumbers=$(find *.mdu | grep -P -o "\d+")
    echo "## Input files start with '$FileName'"
    echo "## Unique identifier is '$Identifier'"
    echo -n "## Found cases: "
    echo $CaseNumbers

    # Loop over all cases
    for Case in $CaseNumbers
    do
        # Wait for a bit
        sleep $(( 10#$Case ))

        # Make use of a queue
        while [ $(jobs -p | wc -l) -ge $Ntasks ]
        do
            sleep 60
        done

        # Prepare
        InputFile="${PWD}/${FileName}${Case}.mdu"
        OutputName="${PWD}/output/${Identifier}_${Case}/FlowFM_map.nc"
        RegridName="${PWD}/output/data_${Identifier}_${Case}.nc"

        # Do computations
        func_computations "$Case" "$InputFile" "$OutputName" "$RegridName" "$Identifier" &
    done

    # Wait for tasks to finish
    wait
    echo "### Finished simulations in $PWD"

    # Return
    cd $BaseDir
done
echo
echo "### Finished all simulations"
echo
