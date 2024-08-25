using ProgressMeter, PyCall, PyPlot, Images, HDF5,NaNStatistics, Statistics, DSP
using _Data
########################################################### parse the arguments
using ArgParse

function string_as_varname_function(s::AbstractString, v::Any)
    s = Symbol(s)
    @eval (($s) = ($v))
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "experiment_filename"
            help = "time stamp, e.g. 20230219_191027"
            required = true
        "server"
            help = "server, e.g. 5"
            required = true
        "experimenter"
            help = "experimenter, e.g. chuyu"
            required = true
        "--analyzer", "-a"
            help = "analyzer, default: chuyu"
            default = "chuyu"
    end

    return parse_args(s)
end


parsed_args = parse_commandline()
# println("Parsed args:")
for (arg,val) in parsed_args
    println("$arg  =>  $val")
    string_as_varname_function(arg, val)
end
########################################################### 
@pyimport numpy

include("../functions/func_map.jl")
include("../functions/func_stat.jl")
include("../functions/func_data.jl")
include("../functions/func_plot.jl")

println("Combine ROI...")


ds_1 = Dataset(experiment_filename, experimenter, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")
ds_save_1 = Dataset(experiment_filename, analyzer, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")

mkpath(data_path(ds_save_1))



try
    NMF_filename = joinpath(data_path(ds_1), "NMF.h5")
    NMF_file = h5open(NMF_filename, "r")
    global z_all = HDF5.readmmap(NMF_file["z_all"])
    global centroid_x_all = HDF5.readmmap(NMF_file["centroid_x_all"])
    global centroid_y_all = HDF5.readmmap(NMF_file["centroid_y_all"])
    global A_all = HDF5.readmmap(NMF_file["A_all"]);
    close(NMF_file)
catch
    NMF_filename = joinpath(data_path(ds_1), "NMF_no_kappa.h5")
    NMF_file = h5open(NMF_filename, "r")
    global z_all = HDF5.readmmap(NMF_file["z_all"])
    global centroid_x_all = HDF5.readmmap(NMF_file["centroid_x_all"])
    global centroid_y_all = HDF5.readmmap(NMF_file["centroid_y_all"])
    global A_all = HDF5.readmmap(NMF_file["A_all"]);
    close(NMF_file)
end

n_rois = size(A_all,2)

X_all = centroid_x_all
Y_all = centroid_y_all
Z_all = z_all;


roi_valid_ratio = [sum(.~isnan.(A_all[:,i])) /prod(size(A_all[:,i])) for i in 1:n_rois];
valid_rois = findall(roi_valid_ratio.>0.3);
component_label = find_same_neuron(X_all[valid_rois], Y_all[valid_rois], Z_all[valid_rois], A_all[:,valid_rois];corr_thres = 0.7);
component_label_small = deal_long_neuron(component_label, Z_all[valid_rois];valid_neurons=valid_rois,A_dF=A_all);


neuron_label = fill(NaN32, n_rois)
neuron_label[valid_rois].=component_label_small;
nr_neuron = Int32.(numpy.nanmax(neuron_label))


h5open(joinpath(data_path(ds_save_1), "NMF_merge.h5"), "w") do file
    file["neuron_label"] = neuron_label
end;

cell_locs = fill(NaN32, nr_neuron, 3)
for which_neuron in 1:nr_neuron
    cell_locs[which_neuron, 1] = mean(centroid_x_all[neuron_label.==which_neuron])
    cell_locs[which_neuron, 2] = mean(centroid_y_all[neuron_label.==which_neuron])
    cell_locs[which_neuron, 3] = mean(z_all[neuron_label.==which_neuron])
end


h5open(joinpath(data_path(ds_save_1), "NMF_merge.h5"), "r+") do file
    file["X_all"] = cell_locs[:, 1]
    file["Y_all"] = cell_locs[:, 2]
    file["Z_all"] = cell_locs[:, 3];
end;


println("Before merging: ", length(valid_rois))
println("After merging: ", numpy.nanmax(neuron_label))


synchrony_offset_sweep, N = get_synchrony_offset_sweep(experiment_filename, server;experimenter=experimenter)

n_sweeps = Int(min(round(N/125), size(A_all,1)-synchrony_offset_sweep))


A_all_merge = fill(NaN32, n_sweeps,nr_neuron)
for which_neuron in 1:nr_neuron
    A_all_merge[:,which_neuron] = numpy.nanmean(A_all[synchrony_offset_sweep+1:synchrony_offset_sweep+n_sweeps,neuron_label.==which_neuron],axis=1)
end

h5open(joinpath(data_path(ds_save_1), "NMF_merge.h5"), "r+") do file
    file["A_all_merge"] = A_all_merge
end;

println("Combine ROI done!")