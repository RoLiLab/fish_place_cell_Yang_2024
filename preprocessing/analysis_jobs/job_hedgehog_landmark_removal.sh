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

./pipeline_neuron_2sessions_rotate.sh 20230409_141947 9 21 20230409_162256 9 24 drew chuyu 0 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230423_103330 8 18 20230423_122640 8 24 jen chuyu 0 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230423_141713 8 22 20230423_155036 8 24 drew chuyu 0 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230423_170228 8 18 20230423_182416 8 23 drew chuyu 0 0 false true $PIPELINE

wait