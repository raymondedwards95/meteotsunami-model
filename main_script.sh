#!/bin/bash
# Runs all simulations in the given folders.
# Also processes the model output.

echo ""
echo "### Running from: main_script.sh"

# Help
TextHelp="
#####
\n
Main script for simulating and visualising.
\n\n
Options:
\n\t -p \t Create parameter files
\n\t -s \t Run simulations, also applies -p
\n\t -a \t Create animations
\n\t -f \t Create figures
\n\t -v \t Create visualisations, applies -f, -v and -t
\n\t -t \t Create figures for theory
\n\n
Folders:
\n\t -r \t Apply options to the reproduction
\n\t -e \t Apply options to the experiments
\n\n
Help:
\n\t -h \t Show this text
\n\n
"

# Maximum number of simultanious tasks
Ntasks=5
Ntasks_half=$(python3 -c "print(f'{max([3, $Ntasks // 3]):0.0f}')")

# Default arguments
run_model=false
create_parameters=false
create_animations=false
create_figures=false
create_theory=false

FolderList=()

# Commandline arguments
echo "## Processing arguments"
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p|--parameters) echo "Got argument '--parameters'" ; create_parameters=true ;;
        -s|--simulations) echo "Got argument '--simulations'" ; run_model=true ;;
        -a|--animations) echo "Got argument '--animations'" ; create_animations=true ;;
        -f|--figures) echo "Got argument '--figures'" ; create_figures=true ;;
        -v|--visualisations) echo "Got argument '--visualisations'" ; create_figures=true ; create_animations=true ; create_theory=true ;;
        -r|--repr) echo "Got argument '--repr'" ; FolderList+=("./reproduction-an-2012") ;;
        -e|--exp) echo "Got argument '--exp'" ; FolderList+=("./thesis") ;;
        -t|--theory) echo "Got argument '--theory'" ; create_theory=true ;;
        -h|--help) echo "Got argument '--help'" ; echo -e $TextHelp ; exit 1 ;;
        *) echo "Unknown parameter $1" ; exit 1 ;;
    esac
    shift
done

echo "## The following sub-tasks will be done:"
if [ "$run_model" = true ] ; then
    create_parameters=true
    echo "# Simulations"
fi
if [ "$create_parameters" = true ] ; then
    echo "# Parameter-files"
fi
if [ "$create_animations" = true ] ; then
    echo "# Animations"
fi
if [ "$create_figures" = true ] ; then
    echo "# Figures"
fi
if [ "$create_theory" = true ] ; then
    echo "# Theory"
fi

# Default list of folders with input files
if [ ${#FolderList[@]} -eq 0 ] ; then
    echo "# Using default folder list"
    FolderList=("./reproduction-an-2012 ./thesis")
fi
echo "## Folders:"
for dir in ${FolderList[@]}
do
    echo "# $dir"
done

sleep 2

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
echo "# Alternative number of simultanious tasks is $Ntasks_half"

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
    local TStart=$(date +%s)

    echo "## Creating parameter files for $LocalIdentifier"
    echo "$(date) - Creating parameter files" >> $LogFile

    # Bathymetry
    echo "# Creating bathymetry files"
    echo "$(date) - Creating bathymetry files" >> $LogFile
    python3 "./create_bathymetry.py" 1> "${LogFolder}/bathymetry_${LocalIdentifier}.log" 2>&1 &

    # Grid
    # NOTE:DOES NOT WORK
    # echo "# Creating grid files"
    # echo "$(date) - Creating grid files" >> $LogFile
    # python3 "./create_grid.py" 1> "${LogFolder}/grid_${LocalIdentifier}.log" 2>&1

    # Observations
    echo "# Creating observation files"
    echo "$(date) - Creating observation files" >> $LogFile
    python3 "./create_observations.py" 1> "${LogFolder}/observation_${LocalIdentifier}.log" 2>&1 &

    # Pressure
    echo "# Creating pressure files"
    echo "$(date) - Creating pressure files" >> $LogFile
    python3 "./create_pressure.py" 1> "${LogFolder}/pressure_${LocalIdentifier}.log" 2>&1 &

    # End
    wait

    local TEnd=$(date +%s)
    local dTTotal=$(python3 -c "print(f'{($TEnd - $TStart) / 60:0.1f}')")

    echo "## Finished creating parameter files in $dTTotal minutes"
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
    echo "# Case $LocalIdentifier case $LocalCase took $dTTotal minutes in total: $dTComputation minutes for D3D and $dTRegrid minutes for regridding"
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

    echo "# Finished animations for case $LocalIdentifier case $LocalCase in $dTTotal minutes"
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

    echo "# Finished figures for case $LocalIdentifier case $LocalCase in $dTTotal minutes"
    echo "$(date) - Finished figures for '$LocalInputFile'" >> $LogFile
}

# Define theory
function func_theory ()
{
    # Start timer
    local TStart=$(date +%s)

    # Create theory
    echo "# Creating figures for theory"
    echo "$(date) - Creating figures for theory" >> $LogFile

    cd ./theory
    echo "# Creating figures in '$PWD'"
    python3 $PWD/theory_part1.py 1> "${LogFolder}/theory_part1.log" 2>&1 &
    python3 $PWD/theory_part2.py 1> "${LogFolder}/theory_part2.log" 2>&1 &
    wait
    cd ..
    echo "# Returned to '$PWD'"

    # End
    local TEnd=$(date +%s)
    local dTTotal=$(python3 -c "print(f'{($TEnd - $TStart) / 60:0.1f}')")

    echo "# Finished figures for theory in $dTTotal minutes"
    echo "$(date) - Finished figures for theory" >> $LogFile
}

# Create theory
if [ "$create_theory" = true ] ; then
    echo "### Create figures theory"
    func_theory
    echo "### Finished theory"
fi

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

    # Create parameter files
    if [ "$create_parameters" = true ] ; then
        echo "### Create parameter files in $PWD"
        func_parameters "$Identifier"
        echo "### Finished parameter files in $PWD"
    fi

    # Do simulations
    if [ "$run_model" = true ] ; then
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
        echo "### Start animations in $PWD"
        for Case in $CaseNumbers
        do
            # Wait for a bit
            sleep $(( 10#$Case ))

            # Make use of a queue
            while [ $(jobs -p | wc -l) -ge $Ntasks_half ]
            do
                sleep 60
            done  # end while loop

            # Create animations
            func_animations "$Case" "$Identifier" &
        done  # end for-loop over cases

        # Wait for simulations to finish
        wait
        echo "### Finished animations in $PWD"
    fi  # end if create_animations

    # Create figures
    if [ "$create_figures" = true ] ; then
        echo "### Start visualising in $PWD"
        for Case in $CaseNumbers
        do
            # Wait for a bit
            sleep $(( 10#$Case ))

            # Make use of a queue
            while [ $(jobs -p | wc -l) -ge 2 ]
            do
                sleep 20
            done  # end while loop

            # Create figures
            func_figures "$Case" "$Identifier" &
        done  # end for-loop over cases

        # Wait for simulations to finish
        wait
        echo "### Finished visualising in $PWD"
    fi  # end if create_figures

    # Return
    cd $BaseDir
done  # end for-loop over folders

echo ""
echo "### Finished all tasks"
echo ""

exit 0
