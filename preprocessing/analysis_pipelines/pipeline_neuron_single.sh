#!/bin/bash

cd ../scripts

DATASET_NAME_1=$1
SERVER_1=$2
GAIN=$3
EXPERIMENTER=$4
ANALYZER=$5
RERUN=$6
WHICH_SCRIPT=$7

if $RERUN = true
then

    ./pipline_registration.sh $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $ANALYZER
    
    julia combine_rois.jl $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER
    
    python baseline_correction_merged.py  $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $GAIN -a $ANALYZER
fi

julia solve_geometry.jl $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER

julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1} --start_index 1801


