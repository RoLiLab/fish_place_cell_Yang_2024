#!/bin/bash

cd ../scripts

DATASET_NAME_1=$1
SERVER_1=$2
GAIN=$3
EXPERIMENTER=$4
ANALYZER=$5
RERUN=$6
WHICH_SCRIPT=$7
DATASET_NAME_2=$8
SERVER_2=$9

if $RERUN = true
then

    ./pipline_registration.sh $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $ANALYZER
    
    julia combine_rois_known.jl $DATASET_NAME_2 $SERVER_2 $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER
    
    python baseline_correction_merged.py  $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $GAIN -a $ANALYZER
fi

julia rotate_experiment_correction.jl $DATASET_NAME_2 $SERVER_2 $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -r 0
julia rotate_experiment_correction.jl $DATASET_NAME_1 $SERVER_1 $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER -r 0

julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1}

julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_2} --skip_shuffle 1


