# fish_place_cell_Yang_2024

This is the repo for paper **A population code for spatial representation in the zebrafish telencephalon**.

Abstract: 

> Spatial learning in teleost fish requires an intact telencephalon, a brain region that contains putative analogs to the mammalian basal ganglia and limbic system (e.g. hippocampus). However, cells fundamental to spatial cognition in mammals (e.g., place cells) have yet to be established in any fish species. In this study, using tracking microscopy to record brain-wide calcium activity in freely swimming larval zebrafish, we compute the spatial information content of each neuron across the brain. Strikingly, in every recorded animal, cells with the highest spatial specificity are enriched in the zebrafish telencephalon. These place cells (PCs) form a population code of space, from which we can decode the animal’s spatial location across time. By continuous recording of population-level activity, we find that the activity manifold of PCs refines and untangles over time. Through systematic manipulation of allothetic and idiothetic cues, we demonstrate that zebrafish PCs integrate multiple sources of information and can flexibly remap to form distinct spatial maps. By analysis of neighborhood distance between PCs across environments, we find evidence for a weakly preconfigured network in the telencephalon. The discovery of zebrafish PCs represents a step forward in our understanding of spatial cognition across species and the functional role of the early vertebrate telencephalon. 

<!-- Cite:

> A population code for spatial representation in the larval zebrafish telencephalon
Chuyu Yang, Lorenz Mammen, Byoungsoo Kim, Meng Li, Drew N Robson, Jennifer M Li
bioRxiv 2023.11.12.566708; doi: https://doi.org/10.1101/2023.11.12.566708 -->

The repo contains:

1. scripts for data preprocessing (`preprocessing`)

2. notebooks for data analysis (`figures`)

3. functions used in data analysis (`functions`)

4. environmental files for dependencies (`env`)

> Some functions in notebooks were imported from `functions` not by using the relative path but absolute path. These need to be adjusted before running the code.


## Dependencies
**Step1: install anaconda for python env management**

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

**Step2: install and activate the saved environment**

```
conda env create -f environment.yml
conda activate place_cell
```

**Step3: install julia packages and change linked python path**

go to `julia` and type `]`, then add necessary packages according to the ones used in notebooks. Packages with underscore (e.g `_Data`) are internal packages from RoLi lab. They are mainly used for file management and are not necessary for reproducing results.

go to `julia`, change the default python path:
```
ENV["PYTHON"] = "$which python" #for example /home/chuyu/miniconda3/envs/place_cell/bin/python
```

then type `]`, and do `build PyCall`


**Step4: make sure cmtk is excutable from the terminal**

https://www.nitrc.org/projects/cmtk

## Corresponding notebooks for each figure

**Figure 1**

a: made in illustrator

b: /figures/lorenz/F2_Decoding/Overview.ipynb "Error by position"

c: /figures/chuyu/figure1/peak_distribution.ipynb

d: /figures/chuyu/figure1/Figure1_maps.ipynb

e: /figures/chuyu/figure1/Figure1_traversal_example.ipynb

f: /figures/chuyu/figure1/Figure1_distribution_cells.ipynb

g: /figures/chuyu/place_cell_overview/place_cell_fraction.ipynb

h: /figures/chuyu/figure1/Figure1_distribution_cells.ipynb

i: /figures/chuyu/figure1/Figure1_maps_edge_center.ipynb
    

**Figure 2**

a-b: /figures/lorenz/F2_Decoding/Overview.ipynb "Error by position"

c-d: /figures/lorenz/F2_Decoding/Overview.ipynb "Example fish"

e: /figures/lorenz/F2_Decoding/Overview.ipynb "error by population"

f: /figures/lorenz/F2_Decoding/Overview.ipynb "error by amount cells"

g: /figures/lorenz/F2_Decoding/Overview.ipynb "minimum cells"

h: /figures/lorenz/F2_Decoding/Overview.ipynb: "overall error distribution"


**Figure 3**

a-b: /figures/lorenz/F3/F3ab_isomap_early_late.ipynb

c,f: /figures/lorenz/F3/F3cf_early_late_specificity.ipynb

d: figures/lorenz/F3/F3d_early_late_decoding.ipynb

e: figures/lorenz/F3/F3e_early_late_field_size.ipynb


**Figure 4**

a, e, i, m: made in illustrator

b-c: /figures/byoungsoo/job_gradient_light_cycle/map_examples_20230111_143729.ipynb

d: /figures/byoungsoo/job_gradient_light_cycle/same_data_tog_new.ipynb

f-g: /figures/chuyu/job_landmark_removal/map_examples_20230419_143319.ipynb

h: /figures/chuyu/job_landmark_removal/same_data_tog.ipynb

j-k: /figures/chuyu/job_boundary_morphing/map_examples_20230224_162111.ipynb

l: /figures/chuyu/job_boundary_morphing/same_data_tog.ipynb

n-o: /figures/chuyu/job_corner_cue_rotation/map_examples_20220406_153842.ipynb

p: /figures/chuyu/job_corner_cue_rotation/same_data_tog_original.ipynb

q: /figures/chuyu/job_corner_cue_rotation/same_data_tog.ipynb


**Figure 5**

a, e, i, m: made in illustrator

b-c: /figures/chuyu/job_hedgehog_rotation_bottom_out/map_examples_20230429_140408.ipynb

d: /figures/chuyu/job_hedgehog_rotation_bottom_out/same_data_tog.ipynb

f-g: /figures/chuyu/job_hedgehog_landmark_removal/map_examples_20230423_103330.ipynb

h: /figures/chuyu/job_hedgehog_landmark_removal/same_data_tog.ipynb

j-k: /figures/chuyu/job_hedgehog_circle/map_examples_20230422_130558.ipynb

l: /figures/chuyu/job_hedgehog_circle/same_data_tog.ipynb

n-o: /figures/chuyu/job_retangular_wall_rotation/map_examples_20230226_122159.ipynb

p: /figures/chuyu/job_retangular_wall_rotation/same_data_tog_original.ipynb

q: /figures/chuyu/job_retangular_wall_rotation/same_data_tog.ipynb



**Figure 6**

a: made in illustrator

b-c: /figures/chuyu/job_corner_cue_circle/map_examples_20220407_152537.ipynb

d: /figures/chuyu/job_corner_cue_circle/same_data_tog.ipynb

e: /figures/chuyu/place_cell_cluster/neighbor_percentage.ipynb


**Supplementary figure 1**

a: /figures/chuyu/place_cell_overview/examples_corner_cue_20220406_111526_different_locs.ipynb

b-c: /figures/chuyu/place_cell_overview/OLS_fit_tel.ipynb

d,f,g,h: /figures/lorenz/EDF1/EDF1dfgh_speed_blurring.ipynb

e: /figures/lorenz/EDF1/EDF1e_firing_variability.ipynb

i-k: /figures/chuyu/multi_modal/examples_multi_modes.ipynb


**Supplementary figure 2**

a: /figures/chuyu/place_cell_overview/geometry_place_cell_percentage.ipynb

b: /figures/chuyu/place_cell_overview/pallium_subpallium_schmatic.ipynb

c: /figures/chuyu/figure1/Figure1_distribution_cells.ipynb

d: /figures/chuyu/place_cell_overview/S3_specificity_cell.ipynb

e: /figures/chuyu/place_cell_overview/S3_spatial_info_cell.ipynb

f: /figures/chuyu/figure1/Figure1_distribution_cells.ipynb

g: /figures/chuyu/place_cell_cluster/example_anatomy.ipynb

h: /figures/chuyu/place_cell_cluster/anatomical_distance_corr.ipynb


**Supplementary figure 3**

a: /figures/chuyu/BVC/BVC_schematics.ipynb

b: made in illustrator

c: /figures/chuyu/BVC/BVC_examples_20230227_165253_with_modelfit.ipynb

d: /figures/chuyu/BVC/BVC_more_examples_calculation.ipynb


**Supplementary figure 4**

a-d: /figures/lorenz/F2_Decoding/Overview.ipynb "Error by position"


**Supplementary figure 5**

a:  /figures/lorenz/F2_Decoding/Overview.ipynb "overall error distribution"

b: /figures/lorenz/F2_Decoding/Overview.ipynb "error by amount cells" with errors_by_amount_240_in_mm_values, errors_by_amount_240_in_mm_keys

c: /figures/lorenz/F2_Decoding/Overview.ipynb "error by long shift"


**Supplementary figure 6**

a: /figures/chuyu/experiments_summary/same_data_tog_with_control_figure_change.ipynb

b-c: /figures/byoungsoo/dropout/place_cell_dropout_merge.ipynb


**Supplementary figure 7**

/figures/chuyu/job_landmark_removal/map_change_distance_cue.ipynb


**Supplementary figure 8**

a-c: /figures/chuyu/experiments_summary/same_data_tog_with_control_figure.ipynb

d: /figures/chuyu/compare_maps/intial_experiments/boundary_morphing/boundary_morph_bestmatch.ipynb

e: /figures/chuyu/experiments_summary/same_data_tog_with_control_figure_change.ipynb

f: /figures/chuyu/job_hedgehog_rotation_bottom_out/clean vs no clean/same_data_tog_clean.ipynb

g: /figures/chuyu/job_hedgehog_rotation_bottom_out/clean vs no clean/same_data_tog_notclean.ipynb


**Supplementary figure 9**

a: made in illustrator

b: /figures/chuyu/job_ABA/example_fish_trajectory.ipynb

c-d: /figures/chuyu/job_ABA/map_examples_AA.ipynb, /figures/chuyu/job_ABA/map_examples_AB.ipynb, /figures/chuyu/job_ABA/map_examples_BA.ipynb

e: /figures/chuyu/job_ABA/same_data_tog_violin.ipynb


**Supplementary figure 10**

/figures/chuyu/experiments_summary/morphing_improvement.ipynb
