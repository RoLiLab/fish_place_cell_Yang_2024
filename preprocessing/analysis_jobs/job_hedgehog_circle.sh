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


./pipeline_neuron_2sessions_change.sh 20230422_160839 8 16 20230422_174803 8 19 lorenz chuyu false $PIPELINE

wait

./pipeline_neuron_2sessions_change.sh 20230422_130558 8 19 20230422_143553 8 23 jen chuyu false $PIPELINE

wait


./pipeline_neuron_2sessions_change.sh 20230617_164214 8 21 20230617_183124 8 24 lorenz chuyu false $PIPELINE

wait