using DelimitedFiles,HDF5, PyCall, PyPlot, NRRD,FileIO, TiffImages,Images, ProgressMeter
using _Data

@pyimport numpy
@pyimport sklearn.decomposition as decomposition
@pyimport matplotlib.patches as patches
@pyimport scipy.ndimage as ndimage

include("../functions/func_data.jl")



experiment_filename = ARGS[1]
server = ARGS[2]
experimenter = ARGS[3]
analyzer = ARGS[4]

ds_save = Dataset(experiment_filename, experimenter, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")
ds_save_analyzer = Dataset(experiment_filename, analyzer, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")

try
    NMF_filename = joinpath(data_path(ds_save), "NMF.h5")
    NMF_file = h5open(NMF_filename, "r")
    try
        global ds_stack_mean = copy(HDF5.readmmap(NMF_file["stack_mean"]))
    catch
        global ds_stack_mean = copy(HDF5.readmmap(NMF_file["ds_stack_mean"]))
    end
    global z_all = HDF5.readmmap(NMF_file["z_all"])
    global centroid_x_all = HDF5.readmmap(NMF_file["centroid_x_all"])
    global centroid_y_all = HDF5.readmmap(NMF_file["centroid_y_all"])

    close(NMF_file)
catch
    NMF_filename = joinpath(data_path(ds_save), "NMF_no_kappa.h5")
    NMF_file = h5open(NMF_filename, "r")
    global z_all = HDF5.readmmap(NMF_file["z_all"])
    global centroid_x_all = HDF5.readmmap(NMF_file["centroid_x_all"])
    global centroid_y_all = HDF5.readmmap(NMF_file["centroid_y_all"])
    # try
    #     global ds_stack_mean = copy(HDF5.readmmap(NMF_file["stack_mean"]))
    # catch
    #     global ds_stack_mean = copy(HDF5.readmmap(NMF_file["ds_stack_mean"]))
    # end
    close(NMF_file)
end

ds_stack_mean[ds_stack_mean.<0] .= 0;

img_save_path = data_path(ds_save_analyzer)
save(joinpath(img_save_path, experiment_filename*"_ds_stack_mean.nrrd"), ds_stack_mean)

println(joinpath(img_save_path, experiment_filename*"_ds_stack_mean.nrrd"))

n_centroid = length(centroid_x_all)
cell_locs = fill(NaN32, n_centroid, 3)

cell_locs[:, 1] = centroid_x_all
cell_locs[:, 2] = centroid_y_all
cell_locs[:, 3] = z_all;

writedlm(joinpath(img_save_path, experiment_filename*"_nmf_loc.txt"), cell_locs)
println(joinpath(img_save_path, experiment_filename*"_nmf_loc.txt"))