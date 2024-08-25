#!/bin/bash

cd ../analysis_pipelines

# DATASET_NAME_1=$1
# SERVER_1=$2
# GAIN_1=$3
# DATASET_NAME_2=$4
# SERVER_2=$5
# GAIN_2=$6
# EXPERIMENTER=$7
# ANALYZER=$8
# RERUN=$9
PIPELINE=place_cell_analysis.jl


./pipeline_neuron_2sessions_change.sh 20220407_152537 4 22 20220407_170908 4 24 jen chuyu false $PIPELINE

wait

./pipeline_neuron_2sessions_change.sh 20220818_123314 7 17 20220818_145358 7 20 chuyu chuyu false $PIPELINE

wait

./pipeline_neuron_2sessions_change.sh 20220818_163902 7 18 20220818_182702 7 21 chuyu chuyu false $PIPELINE

wait

./pipeline_neuron_2sessions_change.sh 20220819_094333 1 17 20220819_114922 1 20 chuyu chuyu false $PIPELINE

wait
