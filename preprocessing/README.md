# Preprocessing for place cell analysis

## Usage
The preprocessing is designed to be used in command-line. To cluster the preprocessing according to the experimental type (number of session and manipulations), the code is organized in a hierachical manner: 

**jobs** -> **pipelines** -> **scipts**

In a **job**, we organize data with the same experimental type together and specifity the details of the data. All the data in one job would use the same pipeline specified by the user, for example:
```
cd ../analysis_pipelines

WHICH_SCRIPT = place_cell_analysis.jl

./pipeline_neuron_single.sh 20220405_171444 4 21 jen chuyu false $WHICH_SCRIPT

wait

./pipeline_neuron_single.sh 20220416_160516 6 22 jen chuyu false $WHICH_SCRIPT

wait

./pipeline_neuron_single.sh 20220417_165530 2 22 jen chuyu false $WHICH_SCRIPT

wait
```



A **pipeline** consists of a series of commands running different scripts. A typical pipleline would follow such a consequence: 
1. combine rois from the same neuron
2. register the average fish brain to the atlas reference brain
3. baseline correction for the flurescent signal
4. solve the chamber geometry
5. calculate place maps and associated qunatities

One example for using the pipeline to process one data that generates the maps for the whole recording session:
```
./pipeline_neuron_single.sh 20231215_123609 6 0 chuyu dogus true place_cell_analysis.jl
```

The passed arguments are documented in `/analysis_pipelines/pipeline_neuron_single.sh`:
```
DATASET_NAME_1=$1
SERVER_1=$2
GAIN=$3
EXPERIMENTER=$4
ANALYZER=$5
RERUN=$6
WHICH_SCRIPT=$7
```

Likewise, another example for using the pipeline to process one data that generates the maps for 3 periods within a recording session:
```
./pipeline_neuron_3periods.sh 20231215_123609 6 0 chuyu dogus 1 3600 3601 7200 7201 10800 false place_cell_analysis.jl
```
The passed arguments are documented in `/analysis_pipelines/pipeline_neuron_3periods.sh`:
```
DATASET_NAME_1=$1
SERVER_1=$2
GAIN=$3
EXPERIMENTER=$4
ANALYZER=$5
SESSION1_START=$6
SESSION1_END=$7
SESSION2_START=$8
SESSION2_END=$9
SESSION3_START=${10}
SESSION3_END=${11}
RERUN=${12}
WHICH_SCRIPT=${13}
```

Other predefined pipelines are in `analysis_pipelines/`.



Details about different **scripts** are listed below:

 `combine_rois.jl`:
combine rois that are close in space, and have similar activities
(generate `NMF_merge.h5`)

 `combine_rois_known.jl`:
similar to combine_rois.jl but with know combining indices, normally used in multi-session experiments to make sure the rois are combined in the same way
(generate `NMF_merge.h5`)

`pipline_registration.sh`:
register to the fish atlas, which consists of several sub-scripts: save_image_coordinate.jl, registration_roli.sh, engert_region.jl

 `save_image_coordinate.jl`:
save the average fish brain images and roi coordinates
(generate `(experiment_filename)_ds_stack_mean.nrrd` and `(experiment_filename)_nmf_loc.txt`)

 `registration_roli.sh`:
register the average fish brain to the atlas reference brain using CMTK
(generate `(experiment_filename)_nmf_loc_(REF_NAME)_registered.txt`)

 `engert_region.jl`:
assign neurons to different brain regions according to atlas masks
(generate `region_roi_bool.h5`)

 `baseline_correction_merged.jl`:
baseline correction for the flurescent signal
(modify `NMF_merge.h5`)

 `solve_geometry.jl`:
solve the chamber geometry from the fish trajectory
(generate `chamber_geometry_(experiment_filename).h5`)

 `solve_geometry_multi.jl`:
similar to solve_geometry.jl but use data from multiple recording sessions
(generate `chamber_geometry_(experiment_filename).h5`)

 `place_cell_analysis.jl`:
calculate place maps and associated qunatities
(generate `for_place_calculation_(loc_roi_filename)_n(n_pos).h5` and `(which_activity)_spatial_info_(start_min)_(end_min)_(loc_roi_filename)_sigma(sigma)_n(n_pos)_A_dF.h5`)

## Data structure
* `NMF_merge.h5`
```
├─ 🔢 A_all_merge #A_all after cell merging
├─ 🔢 A_baseline #estimated baseline for A_all_merge
├─ 🔢 A_dF #baseline subtracted A_all_merge
├─ 🔢 A_dFF #A_dF/(A_baseline+offset), approximates firing rates
├─ 🔢 X_all #X-loc for merged cells
├─ 🔢 Y_all #Y-loc for merged cells
├─ 🔢 Z_all #Z-loc for merged cells
└─ 🔢 neuron_label #labels for merged rois, each number indicates one cell
```

* `chamber_geometry_(experiment_filename).h5`
solved geometry for chamber
```
├─ 🔢 center_loc #center of the chamber
├─ 🔢 chamber_roi #chamber roi
├─ 🔢 countour_array #countour of the chamber
├─ 🔢 heading_r #fish heading after rotation correction (optional)
├─ 🔢 shift #shift for offset correction
├─ 🔢 x_fish #fish x location given the chamber geometry (could be different from the original, for example after rotation correction)
└─ 🔢 y_fish #fish y location given the chamber geometry (could be different from the original, for example after rotation correction)
```

* `for_place_calculation_(loc_roi_filename)_n(n_pos).h5`
binned & digitised fish location
```
├─ 🔢 bool_index #valid indices in place cell analysis (normally after a thresholding on the speed)
├─ 🔢 chamber_roi #chamber roi
├─ 🔢 loc_digital #digitized fish location
├─ 🔢 mask_invalid #mask for the region outside the chamber
├─ 🔢 mask_valid #mask for the region inside the chamber
├─ 🔢 speed_mm_s #fish moving speed
├─ 🔢 valid_moving_indices #indices when the fish is actively moving
├─ 🔢 x_bins
├─ 🔢 x_digital #digitized fish x location (in which x-bins)
├─ 🔢 x_fish_sweep_mean #downsampled fish x location
├─ 🔢 y_bins
├─ 🔢 y_digital #digitized fish y location (in which y-bins)
└─ 🔢 y_fish_sweep_mean #downsampled fish y location
```

* `(which_activity)_spatial_info_(start_min)_(end_min)_(loc_roi_filename)_sigma(sigma)_n(n_pos)_A_dF.h5`
```
├─ 🔢 activity_num_map_all #summed activity for each bin
├─ 🔢 activity_num_map_all_original #summed activity for each bin without gaussian filter
├─ 🔢 bool_index #indices for calculating the maps
├─ 🔢 candidate_place_cell #neurons that pass the stability and specificty threshold
├─ 🔢 dF_mean_all #mean activity for each neuron
├─ 🔢 entropy_all #map entropy
├─ 🔢 entropy_shuffled #map entropy for shuffled data
├─ 🔢 how_many_chunk #number of chunks for stability calculation
├─ 🔢 mask_map_all #number of visits for each bin
├─ 🔢 mask_map_all_original #number of visits for each bin without gaussian filter
├─ 🔢 neuron_valid_percent #percentage of valid values for each neuron
├─ 🔢 place_cell_index #index of place cells
├─ 🔢 place_map_all #place maps
├─ 🔢 place_map_all_original #place maps without gaussian filter
├─ 🔢 place_map_all_thresholded #place maps after thresholding on mask_map_all
├─ 🔢 sigma #ganssian filter variance
├─ 🔢 spatial_info_all #spatial infomration
├─ 🔢 spatial_info_shuffled #spatial infomration for shuffled data
├─ 🔢 specificity
├─ 🔢 specificity_population_z #specificity z-score compared with the population
├─ 🔢 specificity_shuffle_p #specificity percentile compared with shuffled data
├─ 🔢 specificity_shuffle_z #specificity z-score compared with shuffled data
├─ 🔢 stability_all #stability, correaltion among different chunks
└─ 🔢 valid_neurons #which neurons went to the calculation
```


