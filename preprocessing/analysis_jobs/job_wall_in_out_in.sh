#!/bin/bash

cd ../analysis_pipelines

# DATASET_NAME_1=$1
# SERVER_1=$2
# GAIN_1=$3
# DATASET_NAME_2=$4
# SERVER_2=$5
# GAIN_2=$6
# DATASET_NAME_3=$7
# SERVER_3=$8
# GAIN_3=$9
# EXPERIMENTER=${10}
# ANALYZER=${11}
# RERUN=${12}
PIPELINE=place_cell_analysis.jl

./pipeline_neuron_3sessions.sh 20230301_175518 7 15 20230301_191802 7 20 20230301_202336 7 23 chuyu chuyu false $PIPELINE

wait

./pipeline_neuron_3sessions.sh 20230227_155813 1 19 20230227_165253 1 21 20230227_175948 1 24 chuyu chuyu false $PIPELINE

wait

./pipeline_neuron_3sessions.sh 20230301_123207 8 16 20230301_132953 8 19 20230301_142431 8 21 chuyu chuyu false $PIPELINE

wait
