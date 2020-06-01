#!/usr/bin/env bash

#First build and install SICER2
python3 setup.py sdist
source_dist=$(ls dist)
pip3 install ./dist/$source_dist

python ./test/test_scripts/test.py --sicer --recognicer --df --data-dir ./test/test_files --output-dir ./test/test_output --test-dir ./test/expected_output --test-file treatment_1.bed treatment_2.bed --control-file control_1.bed control_2.bed

if [[ $? -eq 0 ]] 
then
	exit 0
else
	exit 1
fi
