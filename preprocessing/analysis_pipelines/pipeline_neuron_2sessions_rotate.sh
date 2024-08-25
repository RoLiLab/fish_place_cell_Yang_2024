#!/bin/bash

cd ../scripts

DATASET_NAME_1=$1
SERVER_1=$2
GAIN_1=$3
DATASET_NAME_2=$4
SERVER_2=$5
GAIN_2=$6
EXPERIMENTER=$7
ANALYZER=$8
ROTATION=$9
SKIP_CORRECTION=${10}
RERUN=${11}
BOTH=${12}
WHICH_SCRIPT=${13}

if $RERUN = true
then
    julia combine_rois.jl $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER
    julia combine_rois_known.jl $DATASET_NAME_1 $SERVER_1 $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER
    
    ./pipline_registration.sh $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $ANALYZER

    python baseline_correction_merged.py  $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $GAIN_1 -a $ANALYZER &
    python baseline_correction_merged.py  $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER $GAIN_2 -a $ANALYZER &  
fi
wait
julia rotate_experiment_correction.jl $DATASET_NAME_1 $SERVER_1 $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER -r $ROTATION -s $SKIP_CORRECTION &
julia rotate_experiment_correction.jl $DATASET_NAME_2 $SERVER_2 $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -r $ROTATION -s $SKIP_CORRECTION &

wait

julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1} --start_index 1801 &
julia $WHICH_SCRIPT $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_2} &

wait

if $BOTH = true
then
    julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_2} --skip_shuffle 1 --start_index 1801 &
    julia $WHICH_SCRIPT $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1} --skip_shuffle 1 &
fi

wait