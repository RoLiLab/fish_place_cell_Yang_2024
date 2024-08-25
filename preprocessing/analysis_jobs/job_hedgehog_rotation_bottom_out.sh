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

./pipeline_neuron_2sessions_rotate.sh 20230406_100925 9 19 20230406_120147 9 24 chuyu chuyu 0 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230407_171016 9 18 20230407_184007 9 23 chuyu chuyu 0 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230429_194403 8 23 20230429_210211 8 24 jen chuyu 0 0 false true $PIPELINE

wait

./pipeline_neuron_2sessions_rotate.sh 20230429_140408 8 16 20230429_153439 8 19 lorenz chuyu 0 0 false true $PIPELINE

wait
