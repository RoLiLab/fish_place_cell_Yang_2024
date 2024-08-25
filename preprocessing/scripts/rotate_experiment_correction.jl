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
        "experimenter"
            help = "experimenter, e.g. chuyu"
            required = true
        "--analyzer", "-a"
            help = "analyzer"
            default = "chuyu"
        "--rotate_angle", "-r"
            help = "rotation angle"
            default = "180"
        "--skip", "-s"
            help = "whether to skip the registration"
            default = "0"       
    end

    return parse_args(s)
end

include("../functions/func_map.jl")
include("../functions/func_stat.jl")
include("../functions/func_data.jl")
include("../functions/func_plot.jl")

println("Rotate experiment correction...")

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

fig_path = joinpath(data_path(ds_save_2), "figures")
mkpath(fig_path)


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


img_bg = img_bg_1;
w = size(img_bg_1, 1)
l = size(img_bg_1, 2);
img_bg_end_1 = img_bg_1[:,:,end]
img_bg_end_2 = img_bg_2[:,:,end];

########################################################### apply rotation and registration
function rotate2d(x, y, θ)
    cosθ = cos(θ)
    sinθ = sin(θ)
    x * cosθ + y * -sinθ, x * sinθ + y * cosθ
end

function rotate2d(x, y, x0, y0, θ)
    x2, y2 = rotate2d(x - x0, y - y0, θ)
    x2 + x0, y2 + y0
end;

if skip == "1"
    rotate_angle = 0
    img_bg_2_inverse = skimage_transform.rotate(img_bg_end_2, rotate_angle);
    shift = [0, 0]
else
    rotate_angle = Int32(floor(parse(Int32, rotate_angle)))
    img_bg_2_inverse = skimage_transform.rotate(img_bg_end_2, rotate_angle);

    shift_calculate, error, diffphase = skimage_registration.phase_cross_correlation(img_bg_2_inverse, img_bg_end_1)
    shift = [shift_calculate[2], shift_calculate[1]]
end


tform = skimage_transform.EuclideanTransform(translation=shift)
img_bg_2_inverse_new = skimage_transform.warp(img_bg_2_inverse, tform);

fig = figure()
imshow(img_bg_end_1', vmax=800)
imshow(img_bg_2_inverse_new', alpha=0.3, cmap= "gray", vmax=600)
fig.savefig(joinpath(fig_path, "$(rotate_angle)_rotation.png"), bbox_inches = "tight",transparent = true,pad_inches = 0)



x_fish_2_r = zeros(Float32, length(x_fish_2))
y_fish_2_r = zeros(Float32, length(y_fish_2))
for i in 1:length(x_fish_2_r)
    x_fish_2_r[i], y_fish_2_r[i] = rotate2d(x_fish_2[i], y_fish_2[i], w/2, l/2, pi*rotate_angle/180)
    x_fish_2_r[i] = x_fish_2_r[i]-shift[2]
    y_fish_2_r[i] = y_fish_2_r[i]-shift[1]
end

heading_2_r = heading_2 .+ pi*rotate_angle/180;
heading_2_r[heading_2_r.>2*pi] .-= 2*pi;

fig = figure()
imshow(img_bg_end_1[:,:,end]', cmap="gray",vmax=400);
scatter(x_fish_2_r[1:125:end],y_fish_2_r[1:125:end], s=1, alpha=0.3)
fig.savefig(joinpath(fig_path, "$(rotate_angle)_rotation_fish_loc.png"), bbox_inches = "tight",transparent = true,pad_inches = 0)


C = vcat(C_1, C_2);
x_fish = vcat(x_fish_1, x_fish_2_r)
y_fish = vcat(y_fish_1, y_fish_2_r);

w = size(img_bg, 1)
l = size(img_bg, 2);
N = length(C)

########################################################### use behavior from two trials to get the chamber roi
countour_array, center_loc,chamber_roi, chamber_roi_x, chamber_roi_y = solve_geometry(x_fish, y_fish, C, img_bg;method="trajectory");


fig = figure(figsize=(15,5))
subplot(1,3,1)
imshow(img_bg[:,:,end], cmap="gray", vmax=500)
subplot(1,3,2)
imshow(chamber_roi)
subplot(1,3,3)
imshow(img_bg[:,:,end], cmap="gray", vmax=500)
imshow(chamber_roi, alpha=0.1, cmap="Reds")
scatter(center_loc[2], center_loc[1])
fig.savefig(joinpath(fig_path, "solved_geometry.png"), bbox_inches = "tight",transparent = true,pad_inches = 0)




h5open(joinpath(data_path(ds_save_1), "chamber_geometry_$(experiment_filename_1).h5"), "w") do file
    file["x_fish"] = x_fish_1
    file["y_fish"] = y_fish_1
    file["shift"] = [0,0]
    file["chamber_roi"] = chamber_roi
    file["center_loc"] = center_loc
    file["countour_array"] = countour_array
    file["heading_r"] = heading_1
end;


h5open(joinpath(data_path(ds_save_2), "chamber_geometry_$(experiment_filename_1).h5"), "w") do file
    file["x_fish"] = x_fish_2_r
    file["y_fish"] = y_fish_2_r
    file["shift"] = shift
    file["chamber_roi"] = chamber_roi
    file["center_loc"] = center_loc
    file["countour_array"] = countour_array
    file["heading_r"] = heading_2_r
end;

println("Rotate experiment correction done!")
