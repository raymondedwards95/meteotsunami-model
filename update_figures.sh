#!/bin/bash

cd ./reproduction-an-2012
bash ./create_figures.sh

cd ..
cd ./thesis
bash ./create_figures.sh

cd ..
cd ./theory
bash ./create_figures.sh

cd ..
exit 0
