#!/usr/bin/env bash

#First build and install SICER2
python3 setup.py sdist
source_dist=$(ls dist)
pip3 install ./dist/$source_dist

sicer_df -t ./test/treatment_1.bed ./test/treatment_2.bed -c ./test/control_1.bed ./test/control_2.bed -s hg38 -o ./test/travisCI/test_output --significant_reads

recognicer_df -t ./test/treatment_1.bed ./test/treatment_2.bed -c ./test/control_1.bed ./test/control_2.bed -s hg38 -o ./test/travisCI/test_output --significant_reads

if python3 ./test/travisCI/compare.py; then
	echo "Test success"
	exit 0
else
	echo "Test failed"
	exit 1
fi
