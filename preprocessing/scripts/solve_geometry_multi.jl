using ProgressMeter, PyCall, PyPlot, Images, HDF5,NaNStatistics, Statistics, DSP, Lasso
using _Data

@pyimport numpy
@pyimport skimage.transform as skimage_transform
@pyimport skimage.registration as skimage_registration 




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
            help = "time stamp of the 1st experiment, e.g. 20230217_160111"
            required = true
        "server_1"
            help = "server of the 1st experiment, e.g. 5"
            required = true
        "experiment_filename_2"
            help = "time stamp of the 2nd experiment, e.g. 20230217_173521"
            required = true
        "server_2"
            help = "server of the 2nd experiment, e.g. 5"
            required = true
        "experiment_filename_3"
            help = "time stamp of the 3rd experiment, e.g. 20230217_173521"
            required = true
        "server_3"
            help = "server of the 3rd experiment, e.g. 5"
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

include("../functions/func_map.jl")
include("../functions/func_stat.jl")
include("../functions/func_data.jl")
include("../functions/func_plot.jl")


parsed_args = parse_commandline()
# println("Parsed args:")
for (arg,val) in parsed_args
    println("$arg  =>  $val")
    string_as_varname_function(arg, val)
end
########################################################### load data
ds_1 = Dataset(experiment_filename_1, experimenter, gethostname() == "roli-$(server_1)" ? "/data" : "/nfs/data$(server_1)")
ds_save_1 = Dataset(experiment_filename_1, analyzer, gethostname() == "roli-$(server_1)" ? "/data" : "/nfs/data$(server_1)")
ds_2 = Dataset(experiment_filename_2, experimenter, gethostname() == "roli-$(server_2)" ? "/data" : "/nfs/data$(server_2)")
ds_save_2 = Dataset(experiment_filename_2, analyzer, gethostname() == "roli-$(server_2)" ? "/data" : "/nfs/data$(server_2)")
ds_3 = Dataset(experiment_filename_3, experimenter, gethostname() == "roli-$(server_3)" ? "/data" : "/nfs/data$(server_3)")
ds_save_3 = Dataset(experiment_filename_3, analyzer, gethostname() == "roli-$(server_3)" ? "/data" : "/nfs/data$(server_3)")


println("Solve geometry (multi session)...")


img_bg_1, x_fish_1, y_fish_1, C_1, heading_1 = h5open(ds_1, "behavior.h5"; raw = true) do file
    read(file, "img_bg"),
    read(file, "fish_anchor_x"),
    read(file, "fish_anchor_y"),
    read(file, "C"),
    read(file, "heading")
end;


img_bg_2, x_fish_2, y_fish_2, C_2, heading_2 = h5open(ds_2, "behavior.h5"; raw = true) do file
    read(file, "img_bg"),
    read(file, "fish_anchor_x"),
    read(file, "fish_anchor_y"),
    read(file, "C"),
    read(file, "heading")
end;


img_bg_3, x_fish_3, y_fish_3, C_3, heading_3 = h5open(ds_3, "behavior.h5"; raw = true) do file
    read(file, "img_bg"),
    read(file, "fish_anchor_x"),
    read(file, "fish_anchor_y"),
    read(file, "C"),
    read(file, "heading")
end;

img_bg = img_bg_1;
w = size(img_bg_1, 1)
l = size(img_bg_1, 2);


C = vcat(C_1, C_2, C_3);
x_fish = vcat(x_fish_1, x_fish_2, x_fish_3)
y_fish = vcat(y_fish_1, y_fish_2, y_fish_3);

w = size(img_bg, 1)
l = size(img_bg, 2);
N = length(C)

########################################################### use behavior from two trials to get the chamber roi
countour_array, center_loc,chamber_roi, chamber_roi_x, chamber_roi_y = solve_geometry(x_fish, y_fish, C, img_bg;method="trajectory");


h5open(joinpath(data_path(ds_save_1), "chamber_geometry_$(experiment_filename_1).h5"), "w") do file
    file["x_fish"] = x_fish_1
    file["y_fish"] = y_fish_1
    file["shift"] = [0,0]
    file["chamber_roi"] = chamber_roi
    file["center_loc"] = center_loc
    file["countour_array"] = countour_array
    file["heading_r"] = heading_1
end;


h5open(joinpath(data_path(ds_save_2), "chamber_geometry_$(experiment_filename_2).h5"), "w") do file
    file["x_fish"] = x_fish_2
    file["y_fish"] = y_fish_2
    file["shift"] = [0,0]
    file["chamber_roi"] = chamber_roi
    file["center_loc"] = center_loc
    file["countour_array"] = countour_array
    file["heading_r"] = heading_2
end;

h5open(joinpath(data_path(ds_save_3), "chamber_geometry_$(experiment_filename_3).h5"), "w") do file
    file["x_fish"] = x_fish_3
    file["y_fish"] = y_fish_3
    file["shift"] = [0,0]
    file["chamber_roi"] = chamber_roi
    file["center_loc"] = center_loc
    file["countour_array"] = countour_array
    file["heading_r"] = heading_3
end;

println("Solve geometry (multi session) done!")