#!/bin/bash

cd ../scripts

DATASET_NAME_1=$1
SERVER_1=$2
GAIN_1=$3
DATASET_NAME_2=$4
SERVER_2=$5
GAIN_2=$6
DATASET_NAME_3=$7
SERVER_3=$8
GAIN_3=$9
EXPERIMENTER=${10}
ANALYZER=${11}
RERUN=${12}
WHICH_SCRIPT=${13}

if $RERUN = true
then
    julia combine_rois.jl $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER
    julia combine_rois_known.jl $DATASET_NAME_1 $SERVER_1 $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER
    julia combine_rois_known.jl $DATASET_NAME_1 $SERVER_1 $DATASET_NAME_3 $SERVER_3 $EXPERIMENTER -a $ANALYZER
    
    ./pipline_registration.sh $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $ANALYZER
    

    python baseline_correction_merged.py  $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER $GAIN_1 -a $ANALYZER &
    python baseline_correction_merged.py  $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER $GAIN_2 -a $ANALYZER &
    python baseline_correction_merged.py  $DATASET_NAME_3 $SERVER_3 $EXPERIMENTER $GAIN_3 -a $ANALYZER &
fi
wait

julia solve_geometry.jl $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER &
julia solve_geometry.jl $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER &
julia solve_geometry.jl $DATASET_NAME_3 $SERVER_3 $EXPERIMENTER -a $ANALYZER &


wait


julia $WHICH_SCRIPT $DATASET_NAME_1 $SERVER_1 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_1} -n 100 --sigma 2 --default_n 0 --start_index 1801
julia $WHICH_SCRIPT $DATASET_NAME_2 $SERVER_2 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_2} -n 100 --sigma 2 --default_n 0
julia $WHICH_SCRIPT $DATASET_NAME_3 $SERVER_3 $EXPERIMENTER -a $ANALYZER -w neuron -l chamber_geometry_${DATASET_NAME_3} -n 100 --sigma 2 --default_n 0