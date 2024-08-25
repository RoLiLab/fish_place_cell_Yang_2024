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


./pipeline_neuron_2sessions_change.sh 20230224_162111 2 17 20230224_171731 2 20 chuyu chuyu false $PIPELINE

wait

./pipeline_neuron_2sessions_change.sh 20230223_135721 9 20 20230223_150240 9 23 chuyu chuyu false $PIPELINE

wait


./pipeline_neuron_2sessions_change.sh 20230618_122231 8 19 20230618_134025 8 22 chuyu chuyu false $PIPELINE

wait

./pipeline_neuron_2sessions_change.sh 20230618_155837 8 19 20230618_170437 8 22 chuyu chuyu false $PIPELINE

wait



