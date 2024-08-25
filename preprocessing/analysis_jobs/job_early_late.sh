#!/bin/bash

cd ../analysis_pipelines

# DATASET_NAME_1=$1
# SERVER_1=$2
# GAIN=$3
# EXPERIMENTER=$4
# ANALYZER=$5
# SESSION1_START=$6
# SESSION1_END=$7
# SESSION2_START=$8
# SESSION2_END=$9
# RERUN=${10}
PIPELINE=place_cell_analysis.jl

./pipeline_neuron_2periods_no_shuffle.sh 20220407_152537 4 0 jen chuyu 1 3600 7201 10800 false $PIPELINE

./pipeline_neuron_2periods_no_shuffle.sh 20220406_111526 9 0 jen chuyu 1 3600 7201 10800 false $PIPELINE

./pipeline_neuron_2periods_no_shuffle.sh 20220407_090156 5 0 jen chuyu 1 3600 7201 10800 false $PIPELINE

./pipeline_neuron_2periods_no_shuffle.sh 20220417_165530 2 0 jen chuyu 1 3600 7201 10800 false $PIPELINE

./pipeline_neuron_2periods_no_shuffle.sh 20220406_153842 9 0 jen chuyu 1 3600 7201 10800 false $PIPELINE

./pipeline_neuron_2periods_no_shuffle.sh 20220405_171444 4 0 jen chuyu 1 3600 7201 10800 false $PIPELINE

./pipeline_neuron_2periods_no_shuffle.sh 20220416_160516 6 0 jen chuyu 1 3600 7201 10800 false $PIPELINE



