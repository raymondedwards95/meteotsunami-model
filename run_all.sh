#!/bin/bash
# Runs all simulations in the given folders.
# Also processes the model output.

# Maximum number of simultanious tasks
Ntasks=5

# Default arguments
run_model=false
create_animations=false
create_figures=false

# Commandline arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--simulations) run_model=true ;;
        -a|--animations) create_animations=true ;;
        -f|--figures) create_figures=true ;;
        -h|--help) echo "Script to run all simulations (option -s), to create all animations (option -a) and to create all figures (option -f)"; exit 1 ;;
        *) echo "Unknown parameter $1"; exit 1 ;;
    esac
    shift
done

echo ""
echo "### Script: run_all.sh"
echo "## Processing arguments"
if [ "$run_model" = true ] ; then
    echo "# Simulations are enabled"
fi
if [ "$create_animations" = true ] ; then
    echo "# Animations are enabled"
fi
if [ "$create_figures" = true ] ; then
    echo "# Figures are enabled"
fi

# List of folders with inputfiles
FolderList=("./reproduction-an-2012 ./thesis")

# Global settings
if [ ! -f ./settings.sh ] ; then
    echo "./settings.sh does not exist. Creating file... Probably it needs to be adjusted manually."
    cp ./setup/settings.sh.template ./settings.sh
    exit 1
fi

BaseDir="${PWD}"
LogFile="${PWD}/last_runs.log"
LogFolder="${PWD}/logs"

mkdir -p "${LogFolder}"

echo ""
echo "### Prepare script"
echo "## Reading './settings.sh'"
source ./settings.sh
echo "# Current directory is '$BaseDir'"
echo "# Logging file is '$LogFile'"
echo "# Path to executables is $DFLOWFM_BIN_PATH"
echo "# Maximum number of simultanious tasks is $Ntasks"

if [ ! -f "$DFLOWFM_BIN_PATH/run_dflowfm.sh" ] ; then
    echo "$DFLOWFM_BIN_PATH/run_dflowfm.sh does not exist. Check ./setup/settings.sh. Exiting..."
    exit 1
fi

# Define parameter creation
function func_parameters ()
{
    # Catch arguments
    local LocalIdentifier=$1

    # Start
    echo "## Creating parameter files for $LocalIdentifier"
    echo "$(date) - Creating parameter files" >> $LogFile

    # Bathymetry
    echo "# Creating bathymetry files"
    echo "$(date) - Creating bathymetry files" >> $LogFile
    python3 "./create_bathymetry.py" 1> "${LogFolder}/bathymetry_${LocalIdentifier}.log" 2>&1

    # Grid
    # NOTE:DOES NOT WORK
    # echo "# Creating grid files"
    # echo "$(date) - Creating grid files" >> $LogFile
    # python3 "./create_grid.py" 1> "${LogFolder}/grid_${LocalIdentifier}.log" 2>&1

    # Observations
    echo "# Creating observation files"
    echo "$(date) - Creating observation files" >> $LogFile
    python3 "./create_observations.py" 1> "${LogFolder}/observation_${LocalIdentifier}.log" 2>&1

    # Pressure
    echo "# Creating pressure files"
    echo "$(date) - Creating pressure files" >> $LogFile
    python3 "./create_pressure.py" 1> "${LogFolder}/pressure_${LocalIdentifier}.log" 2>&1

    # End
    echo "## Finished creating parameter files"
    echo "$(date) - Finished creating parameter files" >> $LogFile
}

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

# Define animations
function func_animations ()
{
    # Catch arguments
    local LocalCase=$1
    local LocalIdentifier=$2

    # Start timer
    local TStart=$(date +%s)

    # Create animations
    echo "# Creating animations for $LocalIdentifier case $LocalCase"
    echo "$(date) - Creating animations for '$LocalInputFile'" >> $LogFile

    python3 $PWD/create_animations.py $LocalCase 1> "${LogFolder}/animations_${LocalIdentifier}_${LocalCase}.log" 2>&1

    # End
    local TEnd=$(date +%s)
    local dTTotal=$(python3 -c "print(f'{($TEnd - $TStart) / 60:0.1f}')")

    echo "# Finished animations for $LocalIdentifier case $LocalCase"
    echo "# Animations for case $LocalIdentifier $LocalCase took $dTTotal minutes"
    echo "$(date) - Finished animations for '$LocalInputFile'" >> $LogFile
}

# Define figures
function func_figures ()
{
    # Catch arguments
    local LocalCase=$1
    local LocalIdentifier=$2

    # Start timer
    local TStart=$(date +%s)

    # Create figures
    echo "# Creating figures for $LocalIdentifier case $LocalCase"
    echo "$(date) - Creating figures for '$LocalInputFile'" >> $LogFile

    python3 $PWD/create_figures.py $LocalCase 1> "${LogFolder}/figures_${LocalIdentifier}_${LocalCase}.log" 2>&1

    # End
    local TEnd=$(date +%s)
    local dTTotal=$(python3 -c "print(f'{($TEnd - $TStart) / 60:0.1f}')")

    echo "# Finished figures for $LocalIdentifier case $LocalCase"
    echo "# Figures for case $LocalIdentifier $LocalCase took $dTTotal minutes"
    echo "$(date) - Finished figures for '$LocalInputFile'" >> $LogFile
}

# Loop over all folders
for Folder in $FolderList
do
    # Move to folder
    cd $Folder
    echo
    echo "### Preparing for tasks in $PWD"

    # Find cases
    FileName=$(find *.mdu | grep -m 1 -P -o "^\D+")
    Identifier=$(find *.mdu | grep -m 1 -P -o "_(\w+)_" | grep -P -o "[^_]*")
    CaseNumbers=$(find *.mdu | grep -P -o "\d+")
    echo "## Input files start with '$FileName'"
    echo "## Unique identifier is '$Identifier'"
    echo -n "## Found cases: "
    echo $CaseNumbers

    # Do simulations
    if [ "$run_model" = true ] ; then
        echo "### Create parameter files in $PWD"
        func_parameters "$Identifier"

        echo "### Start simulations in $PWD"
        for Case in $CaseNumbers
        do
            # Wait for a bit
            sleep $(( 10#$Case ))

            # Make use of a queue
            while [ $(jobs -p | wc -l) -ge $Ntasks ]
            do
                sleep 60
            done  # end while loop

            # Prepare
            InputFile="${PWD}/${FileName}${Case}.mdu"
            OutputName="${PWD}/output/${Identifier}_${Case}/FlowFM_map.nc"
            RegridName="${PWD}/output/data_${Identifier}_${Case}.nc"

            # Do computations
            func_computations "$Case" "$InputFile" "$OutputName" "$RegridName" "$Identifier" &
        done  # end for-loop over cases

        # Wait for simulations to finish
        wait
        echo "### Finished simulations in $PWD"
    fi  # end if run_model

    # Create animations
    if [ "$create_animations" = true ] ; then
        for Case in $CaseNumbers
        do
            # Wait for a bit
            sleep $(( 10#$Case ))

            # Make use of a queue
            while [ $(jobs -p | wc -l) -ge $Ntasks ]
            do
                sleep 60
            done  # end while loop

            func_animations "$Case" "$Identifier" &

        done  # end for-loop over cases

        # Wait for simulations to finish
        wait
        echo "### Finished creating animations in $PWD"
    fi  # end if create_animations

    # Create figures
    if [ "$create_figures" = true ] ; then
        for Case in $CaseNumbers
        do
            # Wait for a bit
            sleep $(( 10#$Case ))

            # Make use of a queue
            while [ $(jobs -p | wc -l) -ge $Ntasks ]
            do
                sleep 60
            done  # end while loop

            func_figures "$Case" "$Identifier" &

        done  # end for-loop over cases

        # Wait for simulations to finish
        wait
        echo "### Finished creating figures in $PWD"
    fi  # end if create figures

    # Return
    cd $BaseDir
done  # end for-loop over folders

echo ""
echo "### Finished all simulations"
echo ""
