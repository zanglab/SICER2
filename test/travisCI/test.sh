#!/usr/bin/env bash

#First build and install SICER2
python3 setup.py sdist
source_dist=$(ls dist)
pip3 install ./dist/$source_dist

python ./test/test_scripts/test.py --sicer --recognicer --df --data_dir ./test/test_files --output_dir ./test/test_output --test_dir ./test/expected_output --test_file treatment_1.bed treatment_2.bed --control_file control_1.bed control_2.bed

if [[ $? -eq 0 ]] 
then
	exit 0
else
	exit 1
fi
