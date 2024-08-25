#!/bin/bash

cd ../scripts

DATASET_NAME_1=$1
SERVER_1=$2
GAIN=$3
EXPERIMENTER=$4
ANALYZER=$5
SESSION1_START=$6
SESSION1_END=$7
SESSION2_START=$8
SESSION2_END=$9
RERUN=${10}
WHICH_SCRIPT=${11}

if $RERUN = true
then
    julia combine_rois.jl $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER
    
    ./pipline_registration.sh $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $ANALYZER

    python baseline_correction_merged.py  $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $GAIN -a $ANALYZER
    
    julia solve_geometry.jl $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER
    
fi


# julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1}

julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1} -s $SESSION1_START -e $SESSION1_END --skip_shuffle 1 --folder place_cell_windows &

julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1} -s $SESSION2_START -e $SESSION2_END --skip_shuffle 1 --folder place_cell_windows &

wait