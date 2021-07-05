#!/bin/bash

## Start
echo "$(date) - Start updating all figures"
echo "$(date) - Start updating all figures" >> last_runs.log

cd ./theory
echo "$(date) --- Creating figures theory"
echo "$(date) --- Creating figures theory" >> ../last_runs.log
bash ./create_figures.sh
echo "$(date) --- Finished figures theory"
echo "$(date) --- Finished figures theory" >> ../last_runs.log


cd ..
cd ./reproduction-an-2012
echo "$(date) --- Creating figures repr"
echo "$(date) --- Creating figures repr" >> ../last_runs.log
bash ./create_figures.sh
echo "$(date) --- Finished figures repr"
echo "$(date) --- Finished figures repr" >> ../last_runs.log

cd ..
cd ./thesis
echo "$(date) --- Creating figures exp"
echo "$(date) --- Creating figures exp" >> ../last_runs.log
bash ./create_figures.sh
echo "$(date) --- Finished figures exp"
echo "$(date) --- Finished figures exp" >> ../last_runs.log

cd ..
echo "$(date) - Finished updating all figures"
echo "$(date) - Finished updating all figures" >> last_runs.log
exit 0
