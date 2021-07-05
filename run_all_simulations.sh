#!/bin/bash

## Start
echo "$(date) - Starting simulations for all cases"
echo "$(date) - Starting simulations for all cases" >> last_runs.log

## Theory
cd ./theory
echo "$(date) --- Creating figures theory"
echo "$(date) --- Creating figures theory" >> ../last_runs.log
bash create_figures.sh

wait
echo "$(date) --- Finished figures theory"
echo "$(date) --- Finished figures theory" >> ../last_runs.log

## Reproduction
cd ..
cd ./reproduction-an-2012
echo "$(date) --- Create input-files repr"
echo "$(date) --- Create input-files repr" >> ../last_runs.log
bash create_all.sh
wait
echo "$(date) --- Finished creating input-files repr"
echo "$(date) --- Finished creating input-files repr" >> ../last_runs.log

echo "$(date) --- Start simulations repr"
echo "$(date) --- Start simulations repr" >> ../last_runs.log
bash run_all.sh
wait
echo "$(date) --- Finished simulations repr"
echo "$(date) --- Finished simulations repr" >> ../last_runs.log

## Thesis
cd ..
cd ./thesis
echo "$(date) --- Create input-files exp"
echo "$(date) --- Create input-files exp" >> ../last_runs.log
bash create_all.sh
wait
echo "$(date) --- Finished creating input-files exp"
echo "$(date) --- Finished creating input-files exp" >> ../last_runs.log

echo "$(date) --- Start simulations exp"
echo "$(date) --- Start simulations exp" >> ../last_runs.log
bash run_all.sh
wait
echo "$(date) --- Finished simulations exp"
echo "$(date) --- Finished simulations exp" >> ../last_runs.log

## End
cd ..
echo "$(date) - Finished simulations for all cases"
echo "$(date) - Finished simulations for all cases" >> last_runs.log
exit 0
