#!/bin/bash

cd ./theory
bash create_figures.sh
wait

cd ..
cd ./reproduction-an-2012
bash create_all.sh
wait
bash run_all.sh
wait

cd ..
cd ./thesis
bash create_all.sh
wait
bash run_all.sh
wait

cd ..
exit 0
