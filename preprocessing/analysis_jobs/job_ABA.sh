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

./pipeline_neuron_3sessions_n_normal.sh 20240128_132056 1 0 20240128_152133 1 0 20240128_164754 1 0 drew chuyu false $PIPELINE

wait

./pipeline_neuron_3sessions_n_normal.sh 20240127_185226 8 0 20240127_202210 8 0 20240127_212910 8 0 jen chuyu true $PIPELINE

wait

