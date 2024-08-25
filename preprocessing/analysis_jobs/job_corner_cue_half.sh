#!/bin/bash

cd ../analysis_pipelines

# DATASET_NAME_1=$1
# SERVER_1=$2
# GAIN=$3
# EXPERIMENTER=$4
# ANALYZER=$5
# RERUN=$6
PIPELINE=place_cell_analysis.jl

./pipeline_neuron_2halves.sh 20220406_111526 9 21 jen chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220407_090156 5 22 jen chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220406_153842 9 21 jen chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220407_152537 4 22 jen chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220417_165530 2 22 jen chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220818_123314 7 17 chuyu chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220818_163902 7 18 chuyu chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220819_094333 1 17 chuyu chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220405_171444 4 21 jen chuyu false $PIPELINE

./pipeline_neuron_2halves.sh 20220416_160516 6 22 jen chuyu false $PIPELINE

