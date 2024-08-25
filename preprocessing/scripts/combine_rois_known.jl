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
        "experiment_filename_1"
            help = "time stamp, e.g. 20230219_191027"
            required = true
        "server_1"
            help = "server, e.g. 5"
            required = true
        "experiment_filename_2"
            help = "time stamp, e.g. 20230219_191027"
            required = true
        "server_2"
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

println("Combine ROI (known)...")


ds_1 = Dataset(experiment_filename_1, experimenter, gethostname() == "roli-$(server_1)" ? "/data" : "/nfs/data$(server_1)")
ds_save_1 = Dataset(experiment_filename_1, analyzer, gethostname() == "roli-$(server_1)" ? "/data" : "/nfs/data$(server_1)")


prev_NMF_filename = joinpath(data_path(ds_save_1), "NMF_merge.h5")
prev_NMF_file = h5open(prev_NMF_filename, "r")
X_all = HDF5.readmmap(prev_NMF_file["X_all"])
Y_all = HDF5.readmmap(prev_NMF_file["Y_all"])
Z_all = HDF5.readmmap(prev_NMF_file["Z_all"])
neuron_label = HDF5.readmmap(prev_NMF_file["neuron_label"])
close(prev_NMF_file)

nr_neuron = Int32.(numpy.nanmax(neuron_label))


ds_2 = Dataset(experiment_filename_2, experimenter, gethostname() == "roli-$(server_2)" ? "/data" : "/nfs/data$(server_2)")
ds_save_2 = Dataset(experiment_filename_2, analyzer, gethostname() == "roli-$(server_2)" ? "/data" : "/nfs/data$(server_2)")

mkpath(data_path(ds_save_2))


h5open(joinpath(data_path(ds_save_2), "NMF_merge.h5"), "w") do file
    file["neuron_label"] = neuron_label
    file["X_all"] = X_all
    file["Y_all"] = Y_all
    file["Z_all"] = Z_all;
end;



try
    NMF_filename = joinpath(data_path(ds_2), "NMF.h5")
    NMF_file = h5open(NMF_filename, "r")
    global A_all = HDF5.readmmap(NMF_file["A_all"]);
    close(NMF_file)
catch
    NMF_filename = joinpath(data_path(ds_2), "NMF_no_kappa.h5")
    NMF_file = h5open(NMF_filename, "r")
    global A_all = HDF5.readmmap(NMF_file["A_all"]);
    close(NMF_file)
end

synchrony_offset_sweep, N = get_synchrony_offset_sweep(experiment_filename_2, server_2;experimenter=experimenter)
n_sweeps = Int(min(round(N/125), size(A_all,1)-synchrony_offset_sweep))

A_all_merge = fill(NaN32, n_sweeps,nr_neuron)
for which_neuron in 1:nr_neuron
    A_all_merge[:,which_neuron] = numpy.nanmean(A_all[synchrony_offset_sweep+1:synchrony_offset_sweep+n_sweeps,neuron_label.==which_neuron],axis=1)
end


h5open(joinpath(data_path(ds_save_2), "NMF_merge.h5"), "r+") do file
    file["A_all_merge"] = A_all_merge
end;

println("Combine ROI (known) done!")