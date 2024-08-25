#!/bin/bash

cd ../analysis_pipelines

# DATASET_NAME_1=$1
# SERVER_1=$2
# GAIN=$3
# EXPERIMENTER=$4
# ANALYZER=$5
# RERUN=$6
PIPELINE=place_cell_analysis_shuffle_all.jl

./pipeline_neuron_single.sh 20220405_171444 4 21 jen chuyu false $PIPELINE

wait

./pipeline_neuron_single.sh 20220416_160516 6 22 jen chuyu false $PIPELINE

wait

./pipeline_neuron_single.sh 20220417_165530 2 22 jen chuyu false $PIPELINE

wait

