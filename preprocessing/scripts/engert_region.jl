using DelimitedFiles,HDF5, PyCall, PyPlot, NRRD,FileIO, TiffImages,Images, ProgressMeter, Statistics,LinearAlgebra
using _Data

@pyimport numpy
# @pyimport sklearn.decomposition as decomposition
# @pyimport matplotlib.patches as patches
@pyimport scipy.ndimage as ndimage

include("../functions/func_map.jl")
include("../functions/func_stat.jl")
include("../functions/func_data.jl")
include("../functions/func_plot.jl")

experiment_filename = ARGS[1]
server = ARGS[2]
experimenter = ARGS[3]
analyzer = ARGS[4]

ds_save = Dataset(experiment_filename, experimenter, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")
ds_save_analyzer = Dataset(experiment_filename, analyzer, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")

cell_locs_original = readdlm(joinpath(data_path(ds_save_analyzer), experiment_filename*"_nmf_loc.txt"));
cell_locs = cell_locs_original[:,1:3];
n_roi = size(cell_locs, 1)

ds_stack_mean = Float32.(load(joinpath(data_path(ds_save_analyzer), experiment_filename*"_ds_stack_mean.nrrd")));

loc_save_path = joinpath(data_path(ds_save_analyzer), "$experiment_filename"*"_engert")
registered_cell_locs_original = readdlm(joinpath(loc_save_path, experiment_filename*"_nmf_loc_engert_registered.txt"));
registered_cell_locs = registered_cell_locs_original[:,1:3];


ds_save_cy_ref = Dataset("20220416_123633", "chuyu", gethostname() == "roli-$(5)" ? "/data" : "/nfs/data$(5)")
# load atlas 
atlas_ref_img = TiffImages.load(data_path(ds_save_cy_ref)*"/Elavl3-H2BRFP.tif");
atlas_ref = convert(Array{Float32}, atlas_ref_img);

atlas_ref_t = fill(NaN32, size(atlas_ref,2), size(atlas_ref,1), size(atlas_ref,3))
for (i, z) in enumerate(1:size(atlas_ref,3))
    atlas_ref_t[:,:,i].= atlas_ref[:,:,z]'
end


atlas_cell_locs = Int32.(round.(registered_cell_locs));
atlas_cell_locs[:,1] .+= 329;
atlas_cell_locs[:,2] .+= 41;

atlas_cell_locs[:,1] .= size(atlas_ref_t,1).-atlas_cell_locs[:,1];
atlas_cell_locs[:,3] .= size(atlas_ref_t,3).-atlas_cell_locs[:,3];

file_path = (gethostname() == "roli-$(7)" ? "/data/jen/data/Z_brain_mask/" : "/nfs/data7/jen/data/Z_brain_mask/")

mask_file = h5open(joinpath(file_path, "Mask.h5"));
region_names = keys(mask_file)
# println(region_names)
# close(mask_file)

chosen_cell_locs = atlas_cell_locs

region_roi_bool = fill(false, n_roi, length(region_names))


println("Get region indices...")
for (i, name) in enumerate(region_names)

    global mask = read(mask_file, name)
    try
        size(mask)
    catch
        mask_key = collect(keys(mask))[1]
        global mask = mask[mask_key]
    end


    whether_in_region = falses(size(chosen_cell_locs,1));
    for i in 1:size(chosen_cell_locs,1)
        x = chosen_cell_locs[i,1]
        y = chosen_cell_locs[i,2]
        z = chosen_cell_locs[i,3]
        if x.<=size(mask,1) && y.<=size(mask,2) && z.<=size(mask,3)&& x>0 && y>0 && z>0
            whether_in_region[i] = mask[x, y, z] #check if it's in any brain region we define
        end
    end
    region_roi_bool[:,i] .= whether_in_region
    
end


h5open(joinpath(data_path(ds_save_analyzer), "region_roi_bool.h5"), "w") do file
    file["region_roi_bool"] = region_roi_bool
    file["region_names"] = region_names
end;