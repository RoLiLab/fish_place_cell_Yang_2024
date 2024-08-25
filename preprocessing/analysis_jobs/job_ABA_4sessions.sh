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

./pipeline_neuron_3sessions_n_normal_0lag.sh 20240127_120909_20240127_132316 1 0 20240127_152329 1 0 20240127_164645 1 0 chuyu chuyu true $PIPELINE

wait

./pipeline_neuron_3sessions_n_normal_0lag.sh 20240128_174842_20240128_182823 8 0 20240128_194250 8 0 20240128_203533 8 0 chuyu chuyu true $PIPELINE

wait
