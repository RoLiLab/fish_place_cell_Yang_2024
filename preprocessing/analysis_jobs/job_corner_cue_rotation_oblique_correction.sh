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
PIPELINE=place_cell_analysis_obilique_correction.jl

./pipeline_neuron_2sessions_rotate.sh 20220406_111526 9 21 20220406_125842 9 24 jen chuyu 180 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20220407_090156 5 22 20220407_104712 5 24 jen chuyu 180 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20220406_153842 9 21 20220406_171558 9 24 jen chuyu 180 0 false true $PIPELINE

wait