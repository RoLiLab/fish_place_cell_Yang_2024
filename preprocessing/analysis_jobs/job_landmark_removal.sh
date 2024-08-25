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
# ROTATION=$9
# SKIP_CORRECTION=${10}
# RERUN=${11}
# BOTH=${12}
PIPELINE=place_cell_analysis.jl


./pipeline_neuron_2sessions_rotate.sh 20230419_143319 9 16 20230419_161007 9 21 chuyu chuyu 0 1 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230420_113755 9 18 20230420_130657 9 21 lorenz chuyu 0 1 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230420_152855 9 14 20230420_165502 9 19 chuyu chuyu 0 1 false true $PIPELINE

wait