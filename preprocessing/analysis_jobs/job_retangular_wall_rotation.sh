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

./pipeline_neuron_2sessions_rotate.sh 20230225_191641 9 21 20230225_210121 9 24 chuyu chuyu 180 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230226_180404 9 21 20230226_190350 9 23 lorenz chuyu 180 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230226_122159 9 21 20230226_141707 9 24 lorenz chuyu 180 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230225_155030 9 22 20230225_174749 9 24 chuyu chuyu 180 0 false true $PIPELINE

wait
