using ProgressMeter, PyCall, PyPlot, Images, HDF5,NaNStatistics, Statistics, DSP, Lasso
using _Data

@pyimport numpy


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

include("../functions/func_map.jl")
include("../functions/func_stat.jl")
include("../functions/func_data.jl")
include("../functions/func_plot.jl")



ds = Dataset(experiment_filename, experimenter, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)") 
ds_save = Dataset(experiment_filename, analyzer, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")


println("Solve geometry...")

C, heading, img_bg, y_fish, x_fish = h5open(ds, "behavior.h5"; raw = true) do file
    read(file, "C"),
    read(file, "heading"), 
    read(file, "img_bg"),
    read(file, "fish_anchor_y"), 
    read(file, "fish_anchor_x")
end;

w = size(img_bg, 1)
l = size(img_bg, 2);
N = length(C)


countour_array, center_loc,chamber_roi, chamber_roi_x, chamber_roi_y = solve_geometry(x_fish, y_fish, C, img_bg;method="trajectory");


center_loc = [(maximum(countour_array[:,1]) + minimum(countour_array[:,1]))/2, (maximum(countour_array[:,2]) + minimum(countour_array[:,2]))/2]


fig_path = joinpath(data_path(ds_save), "figures")
mkpath(fig_path)


fig = figure(figsize=(15,5))
subplot(1,3,1)
imshow(img_bg[:,:,end], cmap="gray", vmax=500)
subplot(1,3,2)
imshow(chamber_roi)
subplot(1,3,3)
imshow(img_bg[:,:,end], cmap="gray", vmax=500)
imshow(chamber_roi, alpha=0.1, cmap="Reds")
scatter(center_loc[2], center_loc[1])
title(experiment_filename)

fig.savefig(joinpath(fig_path, experiment_filename*"_geometry.png"), bbox_inches = "tight",transparent = true,pad_inches = 0)


h5open(joinpath(data_path(ds_save), "chamber_geometry_$(experiment_filename).h5"), "w") do file
    file["chamber_roi"] = chamber_roi
    file["center_loc"] = center_loc
    file["countour_array"] = countour_array
    file["x_fish"] = x_fish
    file["y_fish"] = y_fish
    
end;


println("Solve geometry done!")