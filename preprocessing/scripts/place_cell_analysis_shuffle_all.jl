using ProgressMeter, PyCall, PyPlot, Images, HDF5,NaNStatistics, Statistics, DSP, Lasso, JLD2
using _Data, _Math

@pyimport numpy

include("../functions/func_map.jl")
include("../functions/func_stat.jl")
include("../functions/func_data.jl")
include("../functions/func_plot.jl")

# put all the hyperparameters here
use_gaussian_filter = true
nr_shuffle = 1000
at_least_shift = 1
at_least_frame = 500
at_least_visit = 2



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
            help = "experimenter"
            required = true
        "--analyzer", "-a"
            help = "analyzer"
            default = "chuyu"
        "--which_activity", "-w"
            help = "use which activity to compute the maps, roi or neuron"
            default = "roi"
        "--loc_roi_filename", "-l"
            help = "use which file to compute the maps"
            default = "chamber_geometry"
        "--start_index", "-s"
            help = "start of the map computation"
            default = "-1"
        "--end_index", "-e"
            help = "end of the map computation"
            default = "-1"
        "--start_p"
            help = "start of the map computation (percentage)"
            default = "-1"
        "--end_p"
            help = "end of the map computation (percentage)"
            default = "-1"
        "--n_pos", "-n"
            help = "number of bins"
            default = "90"
        "--sigma"
            help = "sigma for gaussian filter"
            default = "1"
        "--skip_shuffle"
            help = "skip the shuffling part"
            default = "0"
        "--default_n"
            help = "use 60 for small L-bracket, use 90 for large L-bracket"
            default = "1"
        "--exclude_beginning"
            help = "only use the end of the experiment"
            default = "0"
        "--only_end"
            help = "if exclude_beginning, calculate only the last x frames"
            default = "5400"
        "--exclude_end"
            help = "only use the beginning of the experiment"
            default = "0"
        "--only_beginning"
            help = "if exclude_end, calculate only the first x frames"
            default = "5400"
        "--folder"
            help = "the default folder is the data folder"
            default = ""

    end

    return parse_args(s)
end

println("Place cell analysis...")

parsed_args = parse_commandline()
# println("Parsed args:")
for (arg,val) in parsed_args
    println("$arg  =>  $val")
    string_as_varname_function(arg, val)
end

n_pos = parse(Int32, n_pos)
sigma = parse(Int32, sigma)


if loc_roi_filename == "chamber_geometry"
    loc_roi_filename = "chamber_geometry_$(experiment_filename)"
end
########################################################### 

ds = Dataset(experiment_filename, experimenter, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")
ds_save = Dataset(experiment_filename, analyzer, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)")


loc_roi_file = joinpath(data_path(ds_save), loc_roi_filename*".h5")
y_fish_ob, x_fish_ob, chamber_roi = h5open(loc_roi_file) do file
    read(file, "y_fish"), 
    read(file, "x_fish"),
    read(file, "chamber_roi")
end;



rotate_correction = false

if rotate_correction
    orientation_correction_file = h5open(loc_roi_file)
    what_angle = read(orientation_correction_file,"what_angle")
    chamber_roi = read(orientation_correction_file,"chamber_roi_r")
    close(orientation_correction_file)
else 
    what_angle = 0
end


w = size(chamber_roi, 1)
l = size(chamber_roi, 2);
x_fish = zeros(Float32, length(x_fish_ob))
y_fish = zeros(Float32, length(y_fish_ob))
for i in 1:length(x_fish_ob)
    x_fish[i], y_fish[i] = rotate2d(x_fish_ob[i], y_fish_ob[i], w/2, l/2, pi*what_angle/180)
end

N = length(x_fish_ob)

if default_n == "1"
    n_pos = 60
    if w == 5380 # for large L-bracket
        n_pos = 90
    end
end


# Neural data


if which_activity == "roi"
    file = h5open(joinpath(data_path(ds_save), "A_dF_roi.h5"),"r")
    A_dFF = HDF5.readmmap(file["A_dF"]);
    close(file)

    neuron_valid_percent = [sum(.~isnan.(A_dFF[:,i])) /prod(size(A_dFF[:,i])) for i in 1:size(A_dFF,2)];
    valid_neurons = findall(neuron_valid_percent.>=0.3);
    # valid_neurons = 1:size(A_dFF, 2);
elseif which_activity == "neuron"
    file = h5open(joinpath(data_path(ds_save), "NMF_merge.h5"),"r")
    A_dFF = HDF5.readmmap(file["A_dF"]);
    close(file)
    neuron_valid_percent = [sum(.~isnan.(A_dFF[:,i])) /prod(size(A_dFF[:,i])) for i in 1:size(A_dFF,2)];
    valid_neurons = findall(neuron_valid_percent.>=0.3)
end

n_sweeps = size(A_dFF, 1);
n_neurons = size(A_dFF, 2);



if parse(Int32, start_index)<0
    start_index = 1
else
    start_index = maximum([parse(Int32, start_index), 1])
end


if parse(Int32, end_index)<0
    end_index = n_sweeps
else
    end_index = minimum([parse(Int32, end_index), n_sweeps])
end

start_p = parse(Float32, start_p)
end_p = parse(Float32, end_p)

if end_p>0 && end_p<=1
    start_index = maximum([floor(Int32, n_sweeps*start_p), 1])
    end_index = floor(Int32, n_sweeps*end_p)
end    

if exclude_beginning == "1"
    only_end = parse(Int32, only_end)
    start_index = end_index - only_end + 1
    
    @assert start_index > 0
end

if exclude_end == "1"
    only_beginning = parse(Int32, only_beginning)
    end_index = start_index + only_beginning - 1
    end_index = minimum([end_index, n_sweeps])
    # @assert end_index <= n_sweeps
end



start_min = floor(Int32,start_index/120)
end_min = floor(Int32,end_index/120)


save_file_name = "$(which_activity)_spatial_info_$(start_min)_$(end_min)_$(loc_roi_filename)_sigma$(sigma)_n$(n_pos)_A_dF.h5"
save_folder = joinpath(data_path(ds_save), folder)

if end_p>0 && end_p<=1
    save_file_name = "$(which_activity)_spatial_info_$(start_min)_$(end_min)_$(loc_roi_filename)_sigma$(sigma)_n$(n_pos)_A_dF.h5"
    if folder == ""
        folder = "place_cell_windows"
    end
end

if exclude_beginning == "1" || exclude_end == "1"
    save_file_name = "$(which_activity)_spatial_info_$(start_min)_$(end_min)_$(loc_roi_filename)_sigma$(sigma)_n$(n_pos)_A_dF.h5"
end

mkpath(save_folder)

println(joinpath(save_folder, save_file_name))

# Fish location, speed
x_fish_sweep = reshape(x_fish[1:n_sweeps*125], (125, n_sweeps,))
y_fish_sweep = reshape(y_fish[1:n_sweeps*125], (125, n_sweeps,));

x_fish_sweep_mean = dropdims(mean(x_fish_sweep, dims = 1), dims = 1)
y_fish_sweep_mean = dropdims(mean(y_fish_sweep, dims = 1), dims = 1);

x_fish_mm_fit = fit(TrendFilter, x_fish_sweep_mean .* 0.02, 5, 2);
y_fish_mm_fit = fit(TrendFilter, y_fish_sweep_mean .* 0.02, 5, 2);
speed_mm_s = (diff(x_fish_mm_fit.β).^2 .+ diff(y_fish_mm_fit.β).^2).^0.5 ./ 0.5;
speed_mm_s = vcat(speed_mm_s[1], speed_mm_s);
# Moving indices
speed_thr = 0.1
valid_moving_indices = (speed_mm_s .> speed_thr);


# Valid indices for calculation
bool_index = trues(n_sweeps)
if end_index<n_sweeps
    bool_index[end_index+1:end] .= false
end
bool_index[1:start_index-1].=false
bool_index[.!valid_moving_indices].=false;



# Define bins
min_x = 0
min_y = 0
max_x = w
max_y = l

bin_interval = maximum([(max_y-min_y+2)/n_pos,(max_x-min_x+2)/n_pos])

x_bins = collect(min_x-1:bin_interval:min_x+bin_interval*(n_pos)+1);
y_bins = collect(min_y-1:bin_interval:min_y+bin_interval*(n_pos)+1);

x_bins_mid = (x_bins[1:end-1]+x_bins[2:end])/2
y_bins_mid = (y_bins[1:end-1]+y_bins[2:end])/2;

# Digitize data, then we just count the number
x_digital = numpy.digitize(x_fish_sweep_mean, x_bins)
y_digital = numpy.digitize(y_fish_sweep_mean, y_bins);
loc_digital = (y_digital.-1).*n_pos.+x_digital;

chamber_roi_xy = findall(chamber_roi.!=0)
chamber_roi_x = [xy[1] for xy in chamber_roi_xy]
chamber_roi_y = [xy[2] for xy in chamber_roi_xy];
chamber_roi_x_digital = numpy.digitize(chamber_roi_x, x_bins)
chamber_roi_y_digital = numpy.digitize(chamber_roi_y, y_bins);
chamber_roi_digital = (chamber_roi_y_digital.-1).*n_pos.+chamber_roi_x_digital;
mask_valid = MAP.calculate_mask_map_digital(chamber_roi_digital,n_pos).>0;
# mask_valid = erode(dilate(mask_valid))
nr_valid_pos = sum(mask_valid);
mask_invalid = .!mask_valid;
mask_valid_index = findall(mask_valid)
mask_valid_x = [x[1] for x in mask_valid_index]
mask_valid_y = [x[2] for x in mask_valid_index];
mask_valid_x_min = minimum(mask_valid_x)
mask_valid_x_max = maximum(mask_valid_x)
mask_valid_y_min = minimum(mask_valid_y)
mask_valid_y_max = maximum(mask_valid_y);


h5open(joinpath(save_folder, "for_place_calculation_$(loc_roi_filename)_n$(n_pos).h5"), "w") do file
    file["chamber_roi"] = chamber_roi
    file["mask_valid"] = collect(mask_valid)
    file["mask_invalid"] = collect(mask_invalid)
    file["x_bins"] = x_bins
    file["y_bins"] = y_bins
    file["x_digital"] = x_digital
    file["y_digital"] = y_digital
    file["loc_digital"] = loc_digital
    file["x_fish_sweep_mean"] = x_fish_sweep_mean
    file["y_fish_sweep_mean"] = y_fish_sweep_mean
    file["valid_moving_indices"] = collect(valid_moving_indices)
    file["speed_mm_s"] = speed_mm_s
end;


stability_all = fill(NaN32, n_neurons)
how_many_chunk = 2

spatial_info_all = fill(NaN32, n_neurons) #measure for individual neurons
entropy_all = fill(NaN32, n_neurons)
dF_mean_all = fill(NaN32, n_neurons)

activity_num_map_all = fill(NaN32 , n_pos, n_pos, n_neurons)
mask_map_all = fill(NaN32 , n_pos, n_pos, n_neurons)
place_map_all_thresholded = fill(NaN32 , n_pos, n_pos, n_neurons)


activity_num_map_all_original = fill(NaN32 , n_pos, n_pos, n_neurons)
mask_map_all_original = fill(NaN32 , n_pos, n_pos, n_neurons)

println("Calculate maps...")

for neuron_idx in valid_neurons
    neural_activity = A_dFF[:,neuron_idx]
    neural_activity, which_loc = MAP.valid_activity_loc(neural_activity, bool_index,loc_digital)
    if length(neural_activity)<at_least_frame
        continue
    end

    cur_map, mask_map, activity_num_map = MAP.calculate_map_direct(neural_activity, which_loc, n_pos; at_least_visit = at_least_visit, use_gaussian_filter=use_gaussian_filter, sigma=sigma, filter_mask=mask_invalid)
    activity_num_map_all[:,:,neuron_idx] .= activity_num_map
    mask_map_all[:,:,neuron_idx] .= mask_map
    place_map_all_thresholded[:,:,neuron_idx] .= cur_map
    spatial_info_all[neuron_idx] = MAP.spatial_info(cur_map, mask_map)
    entropy_all[neuron_idx] = MAP.map_entropy(cur_map)
    dF_mean_all[neuron_idx] = nanmean(neural_activity)
    
    
    activity_num_map_original = MAP.calculate_activity_map_digital(neural_activity, which_loc,n_pos; use_gaussian_filter=false) #the martix (for each bin) for the summed atcivity
    mask_map_original =MAP.calculate_mask_map_digital(which_loc,n_pos; use_gaussian_filter=false) #the matrix (for each bin) for the number of visits        activity_num_map_all[:,:,neuron_idx] .= activity_num_map
    activity_num_map_all_original[:,:,neuron_idx] .= activity_num_map_original
    mask_map_all_original[:,:,neuron_idx] .= mask_map_original
    
    bool_index_all = generate_chunks(how_many_chunk;n_sweeps=length(neural_activity))
    # bool_index_all = generate_chunks_fine(how_many_chunk;n_sweeps=length(neural_activity))
    map_all = chunk_maps(neural_activity,which_loc, bool_index_all; at_least_visit = at_least_visit, use_gaussian_filter=true, sigma=0.5, filter_mask=mask_invalid)
    stability_all[neuron_idx] = map_stability(map_all)
    
end

place_map_all = activity_num_map_all./mask_map_all;
specificity = spatial_info_all./dF_mean_all;
specificity_population_z = (specificity.-numpy.nanmean(specificity))./numpy.nanstd(specificity);
place_map_all_original = activity_num_map_all_original./mask_map_all_original;

println("Calculate maps done!")

# candidate_place_cell = intersect(findall(stability_all.>0.4), findall(specificity.>0.03))
candidate_place_cell = intersect(findall(specificity_population_z.>=2))

h5open(joinpath(save_folder, save_file_name), "w") do file
    file["bool_index"] = collect(bool_index)
    file["spatial_info_all"] = spatial_info_all
    file["entropy_all"] = entropy_all
    file["place_map_all"] = place_map_all
    file["activity_num_map_all"] = activity_num_map_all
    file["mask_map_all"] = mask_map_all
    file["place_map_all_thresholded"] = place_map_all_thresholded
    file["valid_neurons"] = collect(valid_neurons)
    file["dF_mean_all"] = dF_mean_all
    file["specificity"] = specificity
    file["stability_all"] = stability_all
    file["how_many_chunk"] = how_many_chunk
    file["specificity_population_z"] = specificity_population_z
    file["candidate_place_cell"] = candidate_place_cell
    file["neuron_valid_percent"] = neuron_valid_percent
    file["activity_num_map_all_original"] = activity_num_map_all_original    
    file["mask_map_all_original"] = mask_map_all_original    
    file["place_map_all_original"] = place_map_all_original    
    file["sigma"] = sigma 
end;

if skip_shuffle == "0"

    println("Compare with shuffled...")

    spatial_info_shuffled = fill(NaN32, n_neurons, nr_shuffle) #measure for individual neurons
    entropy_shuffled = fill(NaN32, n_neurons, nr_shuffle)

    shuffle_list = Int32.(numpy.concatenate([-nr_shuffle/2:-1, 1:nr_shuffle/2]));
    for neuron_idx in valid_neurons
    # for neuron_idx in valid_neurons
        neural_activity = A_dFF[:,neuron_idx]
        neural_activity, which_loc = MAP.valid_activity_loc(neural_activity, bool_index,loc_digital)
        if length(neural_activity)<at_least_frame
            continue
        end

        cur_map, mask_map, activity_num_map = MAP.calculate_map_direct(neural_activity, which_loc, n_pos; at_least_visit = at_least_visit, use_gaussian_filter=use_gaussian_filter, sigma=sigma, filter_mask=mask_invalid)

        for i in 1:nr_shuffle
            # randomized_activity = py"shiftshuffle_vector"(neural_activity, at_least_shift)
            perm_idx = numpy.roll(1:length(neural_activity),shuffle_list[i])
            randomized_activity = neural_activity[perm_idx]
            random_activity_num_map  = MAP.calculate_activity_map_digital(randomized_activity, which_loc,n_pos;use_gaussian_filter=use_gaussian_filter, sigma=sigma, filter_mask=mask_invalid)
            random_cur_map= random_activity_num_map./mask_map
            random_cur_map[mask_map .< at_least_visit] .= NaN;
            spatial_info_shuffled[neuron_idx, i] = MAP.spatial_info(random_cur_map, mask_map)
            entropy_shuffled[neuron_idx, i] = MAP.map_entropy(random_cur_map)
        end
    end


    specificity_shuffle_z = (spatial_info_all .- nanmean(spatial_info_shuffled, dims=2)[:,1])./nanstd(spatial_info_shuffled, dims=2)[:,1];
    specificity_shuffle_p = [sum(spatial_info_all[neuron_idx] .> spatial_info_shuffled[neuron_idx,:])/nr_shuffle for neuron_idx in 1:n_neurons];

    place_cell_index = intersect(findall(specificity_population_z.>=3), findall(specificity_shuffle_z.>=5), findall(specificity.>0.01));
    println("Number of place cells: ", length(place_cell_index))

    h5open(joinpath(save_folder, save_file_name), "r+") do file
        file["spatial_info_shuffled"] = spatial_info_shuffled
        file["entropy_shuffled"] = entropy_shuffled
        file["specificity_shuffle_z"] = specificity_shuffle_z    
        file["specificity_shuffle_p"] = specificity_shuffle_p 
        file["place_cell_index"] = place_cell_index
    end;
    
    println("Compare with shuffled done!")
    
end

println("Place cell analysis done!")
