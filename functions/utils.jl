using Distributed
using ProgressMeter


function contrast_stretching(img)
    return (img.-minimum(img))./(maximum(img)-minimum(img))
end

function img_hist(img;  th_l=nothing, th_h=nothing, bins=200)
    img_reshape = reshape(img,size(img,1)*size(img,2));
    if isnothing(th_l)
        th_l = minimum(img_reshape)
    end
    if isnothing(th_h)
        th_h = maximum(img_reshape)
    end  
    img_reshape = img_reshape[findall(img_reshape.<th_h)];
    img_reshape = img_reshape[findall(img_reshape.>th_l)];
    PyPlot.hist(img_reshape,bins);
end

function find_contours(image)
    nbd = 1
    lnbd = 1
    image = Float64.(image)
    contour_list =  Vector{typeof(CartesianIndex[])}()
    done = [false, false, false, false, false, false, false, false]

    # Clockwise Moore neighborhood.
    dir_delta = [CartesianIndex(-1, 0) , CartesianIndex(-1, 1), CartesianIndex(0, 1), CartesianIndex(1, 1), CartesianIndex(1, 0), CartesianIndex(1, -1), CartesianIndex(0, -1), CartesianIndex(-1,-1)]

    height, width = size(image)

    for i=1:height
        lnbd = 1
        for j=1:width
            fji = image[i, j]
            is_outer = (image[i, j] == 1 && (j == 1 || image[i, j-1] == 0)) ## 1 (a)
            is_hole = (image[i, j] >= 1 && (j == width || image[i, j+1] == 0))

            if is_outer || is_hole
                # 2
                border = CartesianIndex[]

                from = CartesianIndex(i, j)

                if is_outer
                    nbd += 1
                    from -= CartesianIndex(0, 1)

                else
                    nbd += 1
                    if fji > 1
                        lnbd = fji
                    end
                    from += CartesianIndex(0, 1)
                end

                p0 = CartesianIndex(i,j)
                detect_move(image, p0, from, nbd, border, done, dir_delta) ## 3
                if isempty(border) ##TODO
                    push!(border, p0)
                    image[p0] = -nbd
                end
                push!(contour_list, border)
            end
            if fji != 0 && fji != 1
                lnbd = abs(fji)
            end

        end
    end

    return contour_list
end

function detect_move(image, p0, p2, nbd, border, done, dir_delta)
    dir = from_to(p0, p2, dir_delta)
    moved = clockwise(dir)
    p1 = CartesianIndex(0, 0)
    while moved != dir ## 3.1
        newp = move(p0, image, moved, dir_delta)
        if newp[1]!=0
            p1 = newp
            break
        end
        moved = clockwise(moved)
    end

    if p1 == CartesianIndex(0, 0)
        return
    end

    p2 = p1 ## 3.2
    p3 = p0 ## 3.2
    done .= false
    while true
        dir = from_to(p3, p2, dir_delta)
        moved = counterclockwise(dir)
        p4 = CartesianIndex(0, 0)
        done .= false
        while true ## 3.3
            p4 = move(p3, image, moved, dir_delta)
            if p4[1] != 0
                break
            end
            done[moved] = true
            moved = counterclockwise(moved)
        end
        push!(border, p3) ## 3.4
        if p3[1] == size(image, 1) || done[3]
            image[p3] = -nbd
        elseif image[p3] == 1
            image[p3] = nbd
        end

        if (p4 == p0 && p3 == p1) ## 3.5
            break
        end
        p2 = p3
        p3 = p4
    end
end

function clockwise(dir)
    return (dir)%8 + 1
end

# rotate direction counterclocwise
function counterclockwise(dir)
    return (dir+6)%8 + 1
end

# move from current pixel to next in given direction
function move(pixel, image, dir, dir_delta)
    newp = pixel + dir_delta[dir]
    height, width = size(image)
    if (0 < newp[1] <= height) &&  (0 < newp[2] <= width)
        if image[newp]!=0
            return newp
        end
    end
    return CartesianIndex(0, 0)
end

# finds direction between two given pixels
function from_to(from, to, dir_delta)
    delta = to-from
    return findall(x->x == delta, dir_delta)[1]
end

function matirx_to_array(matrix)
    dim1 = size(matrix)[1];
    dim2 = size(matrix)[2];
    array = zeros(dim1, dim2);
    for i = 1:dim1
        for j = 1:dim2
            array[i,j] = matrix[i,j];
        end
    end
    return array
end

function distance_from(points, target)
    nr_points = size(points)[1];
    distance = zeros(nr_points);
    
    for index = 1:nr_points
        distance[index] = norm([points[index,1]-target[1], points[index,2]-target[2]]);
    end
    return distance;
end

function angle_from(points, target)
    nr_points = size(points)[1];
    angle_all = Array{Float64, 1}(UndefInitializer(), nr_points);
    
    for index = 1:nr_points
        vector = reshape(points[index,:],1,2).-reshape(target,1,2);
        angle_all[index] = angle(vector[1] + vector[2]*im);
    end
    return angle_all;
end

function angle_between(point, target)
    vector = reshape(point,1,2).-reshape(target,1,2);
    return angle(vector[1] + vector[2]*im);
end

function sliding_window(trace, window_size; mode="variance")
    len_trace = size(trace)[1];
    output_length =  len_trace+1-window_size;
    output = Array{Float64, 1}(UndefInitializer(),output_length);
    if mode == "maximum"
        for i = 1:output_length
            selected= trace[i:window_size+i-1];
            output[i] = maximum(selected);
        end
    elseif mode == "variance"
        for i = 1:output_length
            selected= trace[i:window_size+i-1];
            output[i] = var(selected);
        end
    elseif mode == "sum"
        for i = 1:output_length
            selected= trace[i:window_size+i-1];
            output[i] = sum(selected);
        end
    elseif mode == "dispersal"
        @showprogress for i = 1:output_length
            selected= trace[i:window_size+i-1,:];
            center = mean(selected, dims=1)
            output[i] = sum(distance_from(selected, center));
        end
    end
    return output
end
py"""
import numpy as np
import itertools
import matplotlib.pyplot as plt

def find_consecutive(x):
    return list((g[0][0], len(g)) 
    for key, group in itertools.groupby(enumerate(x), key=lambda v: v[1]) 
    if key 
    for g in (list(group),))

def activity_to_avalanche(summed_activity, min_silence=1):
    avalanche_map=np.zeros_like(summed_activity)
    avalanche_id=1
    if_silence=(summed_activity==0)
    where_silence_period=find_consecutive(if_silence) #find out the silence periods
    where_true_silence_period=[period for period in where_silence_period if period[1]>=min_silence] 
    for which_period in range(1,len(where_true_silence_period)): 
        peirod_now=where_true_silence_period[which_period]
        peirod_past=where_true_silence_period[which_period-1]
        avalanche_map[peirod_past[0]+peirod_past[1]:peirod_now[0]]=avalanche_id 
        avalanche_id+=1

    avalanche_duration=np.zeros(np.max(avalanche_map).astype(int)) 
    avalanche_size=np.zeros(np.max(avalanche_map).astype(int)) 
    avalanche_shape=[]

    for a_index in np.unique(avalanche_map)[1:]:
        a_loc=np.where(avalanche_map==a_index)[0]
        avalanche_shape.append(summed_activity[a_loc])
        avalanche_duration[int(a_index)-1]=len(a_loc)
        avalanche_size[int(a_index)-1]=np.sum(summed_activity[a_loc])
    return avalanche_size,avalanche_duration,avalanche_shape

from math import sqrt

def cercle_circonscrit(T):
    (x1, y1), (x2, y2), (x3, y3) = T
    A = np.array([[x3-x1,y3-y1],[x3-x2,y3-y2]])
    Y = np.array([(x3**2 + y3**2 - x1**2 - y1**2),(x3**2+y3**2 - x2**2-y2**2)])
    if np.linalg.det(A) == 0:
        return False
    Ainv = np.linalg.inv(A)
    X = 0.5*np.dot(Ainv,Y)
    x,y = X[0],X[1]
    r = sqrt((x-x1)**2+(y-y1)**2)
    return (x,y),r

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll

def multicolored_lines(x, y):

    fig, ax = plt.subplots()
    lc = colorline(x, y, ax, cmap='hsv')
    plt.colorbar(lc)
    plt.show()

def colorline(
        x, y, ax, z=None, cmap='copper', norm_min=0, norm_max=1,
        linewidth=1, alpha=1.0):

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(norm_min, norm_max, len(x))

    # Special case if a single number:
    # to check for numerical input -- this is a hack
    if not hasattr(z, "__iter__"):
        z = np.array([z])

    z = np.asarray(z)
    
    segments = make_segments(x, y)
    norm=plt.Normalize(norm_min, norm_max)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    ax.add_collection(lc)

    return lc

def make_segments(x, y):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def shiftshuffle(arr):
    x, y = arr.shape
    rows = np.indices((x,y))[0]
    cols = [np.roll(np.arange(y),np.random.choice(y)) for _ in range(x)]
    return arr[rows, cols]

def shiftshuffle_vector(vec):
    y = len(vec)
    perm_idx = np.roll(np.arange(y),np.random.choice(y))
    return vec[perm_idx]

def crazyshuffle(arr):
    x, y = arr.shape
    rows = np.indices((x,y))[0]
    cols = [np.random.permutation(y) for _ in range(x)]
    return arr[rows, cols]

def crazyshuffle_vector(vec):
    y = len(vec)
    perm_idx = np.random.permutation(y)
    return vec[perm_idx]

def plot_mean_deviation(x0,y0, ax, color='r',plt_raw=True,plt_error=True, nbins=5):
    x=np.array(x0)
    y=np.array(y0)
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    plt_x=(_[1:] + _[:-1])/2
    print('Max of mean is at', round(plt_x[np.argmax(mean)],3),' Max of variance is at', round(plt_x[np.argmax(std)],3))
    if plt_raw==True:
        ax.plot(x, y, 'bo', markersize=1,alpha=0.3, zorder=1, color="gray")
    if plt_error==True:
        errorfill(plt_x, mean, yerr=std, color=color,ax=ax)
    else:
        ax.plot(plt_x, mean, color=color)

def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color, linewidth=3)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)
"""

struct LabelledRanges{T<:Real}
	breakpoints::AbstractArray{T}
	labels::AbstractArray
	isdesc::Bool
end

function LabelledRanges(breakpoints, labels, isdesc=true)
#https://discourse.julialang.org/t/elegant-value-classification-following-a-set-of-rules/14329
    if length(breakpoints) != length(labels)-1
        error("There must be exactly one more label than breakpoints.")
    elseif !issorted(breakpoints)
        error("Breakpoints must be sorted.")
    else
        LabelledRanges(breakpoints, labels, isdesc)
    end
end

getlabel(x::Real, r::LabelledRanges) = r.labels[1+searchsortedlast(r.breakpoints, x, lt=(r.isdesc ? (<=) : (<)))]

function diclists_to_array(diclists)
    all_keys = keys(diclists);
    nr_keys = length(all_keys);
    all_length = [length(diclists[key]) for key in all_keys];
    max_length = maximum(all_length);
    array = zeros(Int64, nr_keys, max_length);
    # array = Array{Any}(nothing,nr_keys, max_length);
    for key in all_keys
        list = each_corner_duration[key];
        array[key,1:length(list)] = list;
    end
    return array
end



function plot_trajectory(fish_position; ax = nothing, fig=nothing,alpha=1,freq=250, plot="plot", s=0.1, t = nothing)
    if isnothing(ax)
        fig, ax = subplots(1,1, figsize=(6,5), dpi=250);
    end
    if isnothing(t)
        t = (1:size(fish_position)[1])/freq/60;
    end
    if plot == "scatter"
        sp = ax.scatter(fish_position[:,1], fish_position[:,2], s=s, c=t,alpha=alpha);
    end
        
    # here t needs to be modified
    if plot == "plot"
        sp = py"colorline"(fish_position[:,1], fish_position[:,2], ax, cmap="viridis", norm_min=minimum(t), norm_max=maximum(t), linewidth=s)
    end
    clb = fig.colorbar(sp,fraction=0.046, pad=0.04)
    clb.ax.set_title("t [min]")
end

function background_mean(img_bg; ifplot=true)
    img_bg_mean = mean(img_bg, dims = 3);
    if length(size(img_bg))>2
        img_bg_mean = dropdims(mean(img_bg, dims = 3), dims = 3);
    else
        img_bg_mean=img_bg;
    end
    if ifplot
        imshow(transpose(img_bg_mean), origin="lower")
    end
    return img_bg_mean
end
;

    
function rotate2d(x, y, θ)
    cosθ = cos(θ)
    sinθ = sin(θ)
    x * cosθ + y * -sinθ, x * sinθ + y * cosθ
end

function rotate2d(x, y, x0, y0, θ)
    x2, y2 = rotate2d(x - x0, y - y0, θ)
    x2 + x0, y2 + y0
end;

function analyse_dispersal(x_pos::Array, y_pos::Array, C::Array, im_freq::Int, time_win::Int, C_lim::Float64)
    "https://github.com/RoLiLab/_BehaviorMetrics.jl.git"
    w = im_freq;
    w2 = round(Int64, w*time_win/2);
    θ = collect(0:0.05:2pi);
    cur_length = zeros(length(θ));
    xy_rotated = zeros(2, 2*w2+1,length(θ));
    temp_trajectory = ones(length(x_pos))*-1;
    
    @showprogress for t = 1:w:length(x_pos)
            idx1 = max(1, t-w2);
            idx2 = min(t+w2, length(x_pos));
            cur_n_frames = idx2-idx1+1;
            if C[t] >= C_lim
                x0, y0 = x_pos[idx1], y_pos[idx1];
                for i = 1:cur_n_frames
                    for j = 1:length(θ)
                            x_r, y_r = rotate2d(x_pos[idx1+i-1], y_pos[idx1+i-1], x0, y0, θ[j])
                            xy_rotated[1,i,j] = x_r
                            xy_rotated[2,i,j] = y_r
                    end 
                end
                for j = 1:length(θ)
                    x1, x2 = extrema(xy_rotated[1,:,j])
                    cur_length[j] = x2-x1
                end
                temp_trajectory[idx1:idx2] .= maximum(cur_length)
            else
               temp_trajectory[idx1:idx2] .= NaN
            end
    end
        
    return(temp_trajectory);
end

function dispersal_2d(x, y)
    center = [mean(x), mean(y)];
    return sum(distance_from(hcat(x, y), center));
end

using PyPlot   # important!
using PyCall
# @pyimport matplotlib.gridspec as gridspec

py"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def imshow_marginal(array, x_min=1, x_max=20):
#https://stackoverflow.com/questions/20525983/matplotlib-imshow-a-2d-array-with-plots-of-its-marginal-densities
    a = np.arange(x_min, x_max+1)
    f = np.arange(1, array.shape[1]+1)
    gs = gridspec.GridSpec(2, 2, width_ratios=[1,3], height_ratios=[3,1])
    ax = plt.subplot(gs[0,1])
    axl = plt.subplot(gs[0,0], sharey=ax)
    axb = plt.subplot(gs[1,1], sharex=ax)

    ax.imshow(array, origin='lower', aspect='auto')
    axl.bar(array.mean(1), a)
    axb.bar(f, array.mean(0))
"""

function dim_one_to_two(one_d, n_row=nothing)
    if isnothing(n_row);
        n_row = maximum(one_d);
    end
    n_column = length(one_d);
    two_d = zeros(n_row, n_column);
    for i in 1:n_column
        if Int(one_d[i]) !=0
            two_d[Int(one_d[i]), i] = 1;
        end
    end
    return two_d
end

function find_silence(speed; plot_state=true,speed_threshold=20,state_threshold_silence = 60*250,state_threshold_active = 10*250,ax=nothing)
    
    whether_active = speed.>speed_threshold
    active_period = py"find_consecutive"(whether_active);
    active_time_duration = active_period[[active_period[i][2]>=state_threshold_active for i =1:length(active_period)]];

    whether_silence=trues(size(whether_active)[1]);

    for i =1:length(active_time_duration)
        whether_silence[active_time_duration[i][1]+1:active_time_duration[i][1]+active_time_duration[i][2]] .= false;
    end


    silence_time_duration = py"find_consecutive"(whether_silence);
    silence_time_duration = silence_time_duration[[silence_time_duration[i][2]>=state_threshold_silence for i =1:length(silence_time_duration)]]

    whether_silence=falses(size(whether_active)[1]);

    for i =1:length(silence_time_duration)
        whether_silence[silence_time_duration[i][1]+1:silence_time_duration[i][1]+silence_time_duration[i][2]] .= true
    end
    
    if plot_state
        if isnothing(ax)
            fig,ax=subplots(1,1,figsize=(15,5))
        end
        ax.plot(whether_silence);

        ax.set_ylabel("Whether silence")
        ax.set_title("silence state threshold = "*string(state_threshold_silence)*", "*"active state threshold = "*string(state_threshold_active)*", "*"speed threshold="*string(speed_threshold))
    end
    return whether_silence, silence_time_duration
end

function string_as_varname_function(s::AbstractString, v::Any)
    s = Symbol(s)
    @eval (($s) = ($v))
end

function h5open_all(file)
    keys_file = keys(file)
    for (index, obj) in enumerate(file)
      data = read(obj)
      string_as_varname_function(keys_file[index], data)
    end
end

function string_as_varname_function(s::AbstractString, v::Any)
    s = Symbol(s)
    @eval (($s) = ($v))
end

function ratio_valid(A)
    return sum(.~isnan.(A))/prod(size(A))
end

# linear interpotor

function sweep_mean(vector::AbstractVector; window_size=ratio_sampling)
    N = length(vector)
    new_N = floor(Int64, N / window_size)
    vector_sweep = reshape(vector[1:new_N*window_size], (window_size, new_N,))
    vector_sweep_mean = dropdims(nanmean(vector_sweep, dims = 1), dims = 1)
    return vector_sweep_mean
end

function sweep_std(vector::AbstractVector; window_size=ratio_sampling)
    N = length(vector)
    new_N = floor(Int64, N / window_size)
    vector_sweep = reshape(vector[1:new_N*window_size], (window_size, new_N,))
    vector_sweep_std = dropdims(nanstd(vector_sweep, dims = 1), dims = 1)
    return vector_sweep_std
end

function filtering(vector::AbstractVector; order=5, λ=2)
    vector_fit = fit(TrendFilter, vector, order, λ);
    return vector_fit.β
end

function speed_2d(x::AbstractVector, y::AbstractVector)
    dx = diff(x)
    dy = diff(y)
    ds_complex = dx + dy.*im
    ds = abs.(ds_complex);
    ds = vcat(ds[1], ds);
    return ds
end

# function clean_nan!(vector::AbstractVector)
#     for i = findall(isnan.(vector))
#         if i==1
#             vector[i]=0
#         else
#             vector[i] = vector[i-1]
#         end
#     end
# end

# function clean_nan(vector::AbstractVector)
#     new_vector = copy(vector)
#     for i = findall(isnan.(vector))
#         if i==1
#             new_vector[i]=0
#         else
#             new_vector[i] = new_vector[i-1]
#         end
#     end
#     return new_vector
# end
py"""
import numpy as np
# ffill along axis 1, as provided in the answer by Divakar
def ffill(arr):
    mask = np.isnan(arr)
    idx = np.where(~mask, np.arange(mask.shape[1]), 0)
    np.maximum.accumulate(idx, axis=1, out=idx)
    out = arr[np.arange(idx.shape[0])[:,None], idx]
    return out

# Simple solution for bfill provided by financial_physician in comment below
def bfill(arr): 
    return ffill(arr[:, ::-1])[:, ::-1]
"""

function clean_nan(matrix)
    matrix_new = py"ffill"(matrix')
    matrix_new = py"bfill"(matrix_new);
    return matrix_new'
end

function distance_to_feature(fish_location::Matrix, fearture_location::Matrix)
    N = size(fish_location)[1]
    fish_x = fish_location[:,1]
    fish_y = fish_location[:,2]
    feature_x = fearture_location[:,1]
    feature_y = fearture_location[:,2]

    dist2feature = zeros(Float32, N);
    nearest_xy = zeros(Float32, 2, N);

    for i = 1:N
        x = fish_x[i]
        y = fish_y[i]
        cur_dist = ((x .- feature_x).^2 + (y .- feature_y).^2).^0.5
        dist2feature[i], min_idx = findmin(cur_dist)
        nearest_xy[1,i] = feature_x[min_idx]
        nearest_xy[2,i] = feature_y[min_idx]
    end
    return dist2feature, nearest_xy
end

function bin_dist_calculate(variable; num_bins=20, bin_size = 1000, sampling_frequency=2)
    half_bin_size = floor(Int64, bin_size)
    variable_min = minimum(variable);
    variable_max = maximum(variable);
    bins = LinRange(variable_min, variable_max, num_bins+1)

    nr_window = length(variable)-bin_size+1;
    array_dist = zeros(num_bins, nr_window);
    mean_hist_all = zeros(nr_window);
    std_hist_all = zeros(nr_window);
    @showprogress for i = 1:nr_window
        chosen_part = variable[i:bin_size+i-1];
        hist_counts = numpy.histogram(chosen_part,bins=bins)[1];

        array_dist[:,i] = hist_counts/sum(hist_counts);
        mean_hist = mean(chosen_part);
        std_hist = std(chosen_part);
        mean_hist_all[i] = mean_hist;
        std_hist_all[i] = std_hist;
    end
    T = half_bin_size:(nr_window+half_bin_size-1)
    T_min = T/sampling_frequency/60
    return array_dist, mean_hist_all, std_hist_all, T_min, bins
end

function bin_dist_plot(variable; ax=nothing, fig=nothing, num_bins=20, bin_size = 1000, sampling_frequency = 2)
    if isnothing(ax)
        fig, ax =subplots(1,1, figsize=(20,5),dpi=250)
    end
    array_dist, mean_hist_all, std_hist_all, T_min, bins = bin_dist_calculate(variable; num_bins=num_bins, bin_size = bin_size, sampling_frequency=sampling_frequency)
    sp = ax.imshow(array_dist, aspect="auto", origin="lower", extent=[minimum(T_min), maximum(T_min),bins[1],bins[end]])

    clb = fig.colorbar(sp,aspect=10)
    clb.ax.set_title("Probability")
    
    yerr = std_hist_all
    ax.plot(T_min,mean_hist_all ,label="mean", color="red",linewidth=3)
    ax.fill_between(T_min, mean_hist_all+std_hist_all, mean_hist_all-std_hist_all, color="orange", alpha=0.25, label="std")
    ax.legend()
end

function load_raw(cur_filename; username="sophie",  host_name="roli-12", which_file=:ir)
    ds = Dataset(cur_filename, username, gethostname() == host_name ? "/data" : "/nfs/data"*host_name[6:end]);
    reader = Reader(ds, which_file);
    N = n_frames(reader);
    print(N)
    return reader
end

function play_video(reader, start, interval, duration; gap=0.01, g=5000)
    for i = start:interval:start+duration
        roi = read(reader, i);
        display(cim(roi./g))
        display(i)
        sleep(gap)
        IJulia.clear_output(true)
    end
end

function play_video(reader, lst; gap=0.01, g=5000)
    for i = lst
        roi = read(reader, i);
        display(cim(roi./g))
        display(i)
        sleep(gap)
        IJulia.clear_output(true)
    end
end


function find_closest_bg(i, img_bg, N)
    nr_bg = size(img_bg, 3)
    bg_interval = floor(Int64, N/size(img_bg, 3))
    bg_timepoints = bg_interval* (1:nr_bg)
    abs_time_dis = abs.(i .- bg_timepoints)
    which_closest = findall(abs_time_dis .== minimum(abs_time_dis))[1]
    img = img_bg[:,:,which_closest]
    return img
end

function calculate_eye_angle(i)
    try
        cur_img = Float32.(read(reader, i));
        img_bg_global = find_closest_bg(i, img_bg)

        cur_img .-= img_bg_global[offset_x[i]+1:offset_x[i]+512, offset_y[i]+1: offset_y[i]+512];
        cur_img[cur_img.< 30] .= 0;
        cur_img[isnan.(cur_img)] .= 0; 
        if _EyeTrackers.track_eye(et, cur_img, Int64(x_fish_fov[i]), Int64(y_fish_fov[i]), Float64(heading[i])) == 0
            return et.theta_left_eye, et.theta_right_eye
        else 
            return NaN,NaN
        end
    catch
        # print(i)
        return NaN,NaN
    end

end


function calculate_eye_angle_diff(i)
    try
        cur_img = Float32.(read(reader, i));
        img_bg_global = find_closest_bg(i, img_bg)

        cur_img .-= img_bg_global[offset_x[i]+1:offset_x[i]+512, offset_y[i]+1: offset_y[i]+512];
        cur_img .-= img_bg_local3;
        cur_img[cur_img.< 30] .= 0;
        cur_img[isnan.(cur_img)] .= 0; 
        if _EyeTrackers.track_eye(et, cur_img, Int64(x_fish_fov[i]), Int64(y_fish_fov[i]), Float64(heading[i])) == 0
            return et.theta_left_eye, et.theta_right_eye
        else 
            return NaN,NaN
        end
    catch
        # print(i)
        return NaN,NaN
    end

end

py"""
import numpy as np
from numba import njit, prange
from tqdm import tqdm
@njit(parallel=True)
def Normalization(fluoTrace,tauDecay,fps, k):
    numFrames=len(fluoTrace)
    twdw=max([15,k*tauDecay])
    wdw=int(round(fps*twdw))
    smoothBaseline=np.zeros_like(fluoTrace)

    temp=np.zeros(numFrames-2*wdw)
    for i in range(wdw,numFrames-wdw):
        temp[i-wdw]=np.nanpercentile(fluoTrace[i-wdw:i+wdw],8)
    smoothBaseline[:wdw]=temp[0]*np.ones(wdw)
    smoothBaseline[wdw:numFrames-wdw]=temp
    smoothBaseline[numFrames-wdw:]=temp[-1]*np.ones(wdw)

    deltaFoF=(fluoTrace-smoothBaseline)/smoothBaseline
    return deltaFoF

@njit(parallel=True)
def window_baseline(fluoTrace, wdw, p=10):
    numFrames=len(fluoTrace)
    smoothBaseline=np.zeros_like(fluoTrace)
    temp=np.zeros(numFrames-2*wdw)
    for i in range(wdw,numFrames-wdw):
        temp[i-wdw]=np.nanpercentile(fluoTrace[i-wdw:i+wdw],p)
    smoothBaseline[:wdw]=temp[0]*np.ones(wdw)
    smoothBaseline[wdw:numFrames-wdw]=temp
    smoothBaseline[numFrames-wdw:]=temp[-1]*np.ones(wdw)
    return smoothBaseline
"""

function show_location(which_neuron; s= 1)
    test_z_idx = z_all[which_neuron];
    scatter(centroid_x_all, centroid_y_all,s=1, color ="gray")
    # img_mean_all = nanmean(stack_mean[:, :, unique(test_z_idx)], dim = 3);
    # imshow(img_mean_all', "bone")
    scatter(centroid_x_all[which_neuron], centroid_y_all[which_neuron],s=s, color ="r", alpha = 0.2)

end

function show_location(which_neuron, c; clb_name = "", cmap="bwr")
    test_z_idx = z_all[which_neuron];
    scatter(centroid_x_all, centroid_y_all,s=1, color ="gray")
    # img_mean_all = nanmean(stack_mean[:, :, unique(test_z_idx)], dim = 3);
    # imshow(img_mean_all', "bone")
    scatter(centroid_x_all[which_neuron], centroid_y_all[which_neuron], s=1, c = c, cmap=cmap)
    clb = colorbar()
    clb.ax.set_title(clb_name)

end

function show_location(which_neuron, ax::PyObject ; s= 1)
    test_z_idx = z_all[which_neuron];
    ax.scatter(centroid_x_all, centroid_y_all,s=1, color ="gray")
    # img_mean_all = nanmean(stack_mean[:, :, unique(test_z_idx)], dim = 3);
    # imshow(img_mean_all', "bone")
    ax.scatter(centroid_x_all[which_neuron], centroid_y_all[which_neuron],s=s, color ="r", alpha = 0.2)

end

function show_location(which_neuron, c, ax::PyObject , fig; clb_name = "", cmap="bwr")
    test_z_idx = z_all[which_neuron];
    ax.scatter(centroid_x_all, centroid_y_all,s=1, color ="gray")
    # img_mean_all = nanmean(stack_mean[:, :, unique(test_z_idx)], dim = 3);
    # imshow(img_mean_all', "bone")
    ax.scatter(centroid_x_all[which_neuron], centroid_y_all[which_neuron], s=1, c = c, cmap=cmap)
    clb = fig.colorbar()
    clb.ax.set_title(clb_name)

end
function generate_map(neural_activity, x, y, x_bins, y_bins, valid_index)
    H_weighted, xedges, yedges = numpy.histogram2d(x[valid_index], y[valid_index], [x_bins, y_bins], weights = neural_activity[valid_index]);
    H, xedges, yedges = numpy.histogram2d(x[valid_index], y[valid_index], [x_bins, y_bins]);
    return H_weighted./H;   
end

function generate_map(neural_activity, x, x_bins,valid_index)
    H_weighted, xedges, yedges = numpy.histogram(x[valid_index], x_bins, weights = neural_activity);
    H, xedges, yedges = numpy.histogram(x[valid_index], x_bins);
    return H_weighted./H;   
end

function normalizing_map(map)
    min_map = numpy.nanmin(map)
    normalized_map = map .- min_map
    sum_map = sum(normalized_map)
    if sum_map != 0
        normalized_map = normalized_map./sum_map
        return normalized_map
    else 
        return []
    end
end

function event_trigger(events_index, response, interval; ax=nothing)
    if isnothing(ax)
        fig, ax = subplots(figsize=(10,5))
    end
    event_response_all = []
    for event in events_index
        ax.plot(-interval:interval, response[event-interval:event+interval], color="gray", alpha=0.2)
        append!(event_response_all,[response[event-interval:event+interval]])
    end
    averaged_response, _ = mean_std(event_response_all)
    ax.plot(-interval:interval, averaged_response, color="black")
end
    
using NaNStatistics

function list_array_to_2d(list_array)
    nr_array = size(list_array)[1]
    length_each = [size(array)[1] for array in list_array]
    new_array = fill(NaN, (nr_array,maximum(length_each)))
    for i = 1:nr_array
        new_array[i,1:length_each[i]] = list_array[i]
    end
    return new_array
end

function mean_std(list_array)
    new_array = list_array_to_2d(list_array)
    mean_array = nanmean(new_array,dims=1)[:]
    std_array = nanmean(new_array, dims=1)[:]
    return mean_array, std_array
end

function hist_m(m, bins=100)
    hist(reshape(m ,length(m)), bins=bins)
end

@pyimport numpy
function correlate_1d(trace1, trace2=nothing; shift = 50, exclude_dx_0 = true)
    """
    trace shape: T*nr_variable
    """
    if isnothing(trace2)
        trace2 = copy(trace1)
    end
    #initialize the output, reformat the map
    len_trace1 = size(trace1, 1)
    len_trace2 = size(trace2, 1)
    auto_correlate = fill(NaN32, 2*shift+1)
    trace1_new = fill(NaN32, len_trace1+2*shift)
    trace1_new[shift+1:len_trace1+shift] .=trace1;
    for (i, dx) in enumerate(-shift:shift)
        trace2_new = fill(NaN32, len_trace2+2*shift)
        
        trace2_new[shift+1+dx:len_trace2+shift+dx] .=trace2
        valid_index = findall((.!isnan.(trace1_new)).*(.!isnan.(trace2_new))) #find the overlapping pixel indices
        if length(valid_index)>0
            auto_correlate[i] = numpy.corrcoef(trace1_new[valid_index],trace2_new[valid_index])[1,2] #calculate correlation
        end
    end
    if exclude_dx_0
        auto_correlate[shift+1] = NaN32
    end
    return auto_correlate
end


function correlate_nd(trace1, trace2=nothing, shift = 50, exclude_dx_0 = true)
    """
    trace shape: T*nr_variable
    under construction
    """
    if isnothing(trace2)
        trace2 = copy(trace1)
    end
    #initialize the output, reformat the map
    len_trace1 = size(trace1, 1)
    len_trace2 = size(trace2, 1)
    auto_correlate = fill(NaN32, 2*shift+1)
    trace1_new = fill(NaN32, len_trace1+2*shift)
    trace1_new[shift+1:len_trace1+shift] .=trace1;
    for (i, dx) in enumerate(-shift:shift)
        trace2_new = fill(NaN32, len_trace2+2*shift)
        
        trace2_new[shift+1+dx:len_trace2+shift+dx] .=trace2
        valid_index = findall((.!isnan.(trace1_new)).*(.!isnan.(trace2_new))) #find the overlapping pixel indices
        if length(valid_index)>0
            auto_correlate[i] = numpy.corrcoef(trace1_new[valid_index],trace2_new[valid_index])[1,2] #calculate correlation
        end
    end
    if exclude_dx_0
        auto_correlate[shift+1] = NaN32
    end
    return auto_correlate
end

function whether_in(vector, collection)
    return [x in collection for x in vector]
end

py"""
def smooth(a,WSZ):
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),"valid")/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))
"""

function PCA_ICA(A_dF, k=2, tol=0.0001)
    """
    A_dF: n_sweeps * n_neurons
    k: ratio for the elbow cutoff
    """
    N = size(A_dF, 1)
    P = numpy.percentile(reshape(A_dF, length(A_dF)),90);
    A_normalized = A_dF.*(10000/P);
    U1, s1, Vh1 = extmath.randomized_svd(A_normalized, n_components=minimum([1000,N]));
    Spectrum=log.(s1);
    Diffvalue=-(maximum(Spectrum)-minimum(Spectrum))/1000;
    DiffFunc=diff(py"smooth"(Spectrum,23));
    NPCfilt = minimum([findfirst(DiffFunc.>Diffvalue), N])

    uu = U1[:, 1:NPCfilt]
    ss = numpy.diag(s1[1:NPCfilt])
    vv = Vh1[1:NPCfilt,:];
    vv[abs.(vv).<0.4*std(vv)].=0;

    stddevs = std(A_normalized-uu*ss*vv, dims=1)
    stddevs[stddevs.<1].=1;

    A_normalized_s = A_normalized./repeat(stddevs, N, 1);
    U2, s2, Vh2 = extmath.randomized_svd(A_normalized_s, n_components=minimum([1000,N]));

    Spectrum=log.(s2);
    Diffvalue=-(maximum(Spectrum)-minimum(Spectrum))/1000;
    DiffFunc=diff(py"smooth"(Spectrum,23));
    Npc = minimum([k*findfirst(DiffFunc.>Diffvalue),N])

    U3 = U2[:, 1:Npc]
    s3 = s2[1:Npc]
    Vh3 = Vh2[1:Npc,:]
    pc_maps = Vh3

    transformer = decomposition.FastICA(n_components=Npc, tol=tol)
    icasig = transformer.fit_transform(pc_maps');
    A = transformer.mixing_;

    GM=icasig;
    std_GM = std(GM, dims=1);
    GMv = GM./repeat(std_GM, size(GM, 1), 1)

    GMzp=GMv.-2;
    GMzp[GMzp.<0].=0;
    GMzpm=mean(GMzp, dims=1)

    GMzn=GMv.+2;
    GMzn[GMzp.>0].=0;
    GMznm=mean(GMzn, dims=1)

    GMs=GM.*repeat(sign.(GMzpm+GMznm),size(GM,1),1);

    TS=U3*numpy.diag(s3)*A;
    TSs=TS.*repeat(sign.(GMzpm+GMznm),size(TS,1),1);

    var_TS = var(TS, dims=1)[1,:]
    p = sortperm(var_TS)
    p = reverse(p)
    TSo = TSs[:,p]
    ic_maps = GMs[:,p];
    return pc_maps, ic_maps, TSo
end

py"""
def smooth(a,WSZ):
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),"valid")/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))
"""

function corr_nan(trace1, trace2; at_least_overlap = 2)
    correlation = NaN
    valid_index = findall((.!isnan.(trace1)).*(.!isnan.(trace2))) #find the overlapping pixel indices
    if length(valid_index) >= at_least_overlap
        correlation = numpy.corrcoef(trace1[valid_index],trace2[valid_index])[1,2] #calculate correlation
    end
    return correlation
end


function corr_nan_matrix(matrix; at_least_overlap = 2)
    nr_variable = size(matrix, 2)
    corr_matrix = numpy.eye(nr_variable)
    i_j_combo = [(i,j) for i in 1:nr_variable for j in i+1:nr_variable]
    for (i,j) in i_j_combo
        vector_i = matrix[:, i]
        vector_j = matrix[:, j]
        corr_value = corr_nan(vector_i, vector_j; at_least_overlap = 2)
        corr_matrix[i,j] = corr_value
        corr_matrix[j,i] = corr_value
    end
    return corr_matrix
end

function correct_angle(angle_vec)
    angle_vec_new = copy(angle_vec)
    larger_than_pi = (angle_vec.>=pi);
    smaller_than_mpi = (angle_vec.<=-pi);
    angle_vec_new[larger_than_pi] = -(2*pi .- angle_vec[larger_than_pi]);
    angle_vec_new[smaller_than_mpi] = 2*pi .+ angle_vec[smaller_than_mpi];
    return angle_vec_new
end





function cell_footprint(which_neuron)
    W = 32
    H=32
    
    img_width = size(ds_stack_mean, 1)
    img_height = size(ds_stack_mean, 2)

    roi_x1 = roi_x_idx_all[which_neuron]
    roi_x2 = min(img_width, roi_x1+W-1)
    roi_x2 == img_width ? roi_x1 = img_width-W+1 : roi_x1 = roi_x1 
    roi_y1 = roi_y_idx_all[which_neuron]
    roi_y2 = min(img_height, roi_y1+H-1)
    roi_y2 == img_height ? roi_y1 = img_height-H+1 : roi_y1 = roi_y1 

    cur_img_x = round(Int64, roi_centroid_x_all[which_neuron] + roi_x1 -1)
    cur_img_y = round(Int64,roi_centroid_y_all[which_neuron] + roi_y1 -1)
    footprint = S_all[:, which_neuron]
    # footprint ./= maximum(footprint);
    cell_loc_map = fill(NaN32, img_width, img_height)
    cell_loc_map[roi_x1:roi_x2, roi_y1:roi_y2].= reshape(footprint, W, H)
    return cell_loc_map
end

function smooth_variable(variable, window)
    smooth_variable = zeros(Float32, length(variable))
    for i in eachindex(variable)
        st = max(1,i-window)
        ed = min(i+window,length(variable))
        smooth_variable[i] = mean(variable[st:ed])
    end
    return smooth_variable
end;



# py"""
# import numpy as np
# def smooth_scatter(x,y,bins =10):
#     n, _ = np.histogram(x, bins=bins)
#     sy, _ = np.histogram(x, bins=bins, weights=y)
#     sy2, _ = np.histogram(x, bins=bins, weights=y*y)
#     mean = sy / n
#     regress_y = mean
#     regress_x=(_[1:] + _[:-1])/2
#     std = np.sqrt(sy2/n - mean*mean)
#     return regress_x, regress_y, std
# """


# py"""
# import numpy as np
# def smooth_scatter_sem(x,y,bins =10):
#     n, _ = np.histogram(x, bins=bins)
#     sy, _ = np.histogram(x, bins=bins, weights=y)
#     sy2, _ = np.histogram(x, bins=bins, weights=y*y)
#     mean = sy / n
#     regress_y = mean
#     regress_x=(_[1:] + _[:-1])/2
#     sem = np.sqrt(sy2/n - mean*mean)/np.sqrt(n)
#     return regress_x, regress_y, sem
# """


py"""
import numpy as np
import scipy.stats as stats
def smooth_scatter_sem(x,y,bins =10):
    n, _, binnumber = stats.binned_statistic(x, y, bins =bins, statistic='count')
    mean, _, binnumber = stats.binned_statistic(x, y, bins =bins, statistic='mean')
    std, _, binnumber = stats.binned_statistic(x, y, bins =bins, statistic='std')
    regress_y = mean
    regress_x=(_[1:] + _[:-1])/2
    return regress_x, regress_y, std/np.sqrt(n), n
"""

py"""
import numpy as np
import scipy.stats as stats
def smooth_scatter(x,y,bins =10):
    n, _, binnumber = stats.binned_statistic(x, y, bins =bins, statistic='count')
    mean, _, binnumber = stats.binned_statistic(x, y, bins =bins, statistic='mean')
    std, _, binnumber = stats.binned_statistic(x, y, bins =bins, statistic='std')
    regress_y = mean
    regress_x=(_[1:] + _[:-1])/2
    return regress_x, regress_y, std, n
"""



