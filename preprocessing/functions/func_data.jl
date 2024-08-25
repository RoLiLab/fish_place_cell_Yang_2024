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

function h5map_all(file)
    keys_file = keys(file)
    for (index, obj) in enumerate(file)
       data = HDF5.readmmap(file[keys_file[index]])
      string_as_varname_function(keys_file[index], data)
    end
end

function clean_nan(vector::AbstractVector)
    new_vector = copy(vector)
    for i = findall(isnan.(vector))
        if i==1
            new_vector[i]=0
        else
            new_vector[i] = new_vector[i-1]
        end
    end
    return new_vector
end


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


function get_synchrony_offset_sweep(experiment_filename,server; experimenter="jen")
    ds = Dataset(experiment_filename, experimenter, gethostname() == "roli-$(server)" ? "/data" : "/nfs/data$(server)") 

    input_filename = joinpath(data_path(ds), "registration_experiment.h5")
    p_rigid = h5read(input_filename, "p_rigid");
    c_rigid = h5read(input_filename, "c_rigid");
    
    heading = h5open(ds, "behavior.h5", "r"; raw = true) do file
        read(file, "heading")
    end;
    fl_heading = mod.(p_rigid[4, :] .+ 4 * pi .- pi /2, 2 * pi);
    fl_heading[c_rigid .< 0.6] .= NaN; 

    for i = 1:length(fl_heading)
        if isnan(fl_heading[i])
            fl_heading[i] = i == 1 ? fl_heading[findfirst(isfinite.(fl_heading))] : fl_heading[i - 1]
        end
    end

    fl_heading = fl_heading[1:2:end]; 
    nir_heading = heading[1:5:end];
    cur_n = min(length(fl_heading), length(nir_heading))
    fl_heading2 = fl_heading[1:cur_n] .- mean(fl_heading[1:cur_n])
    fl_heading2 ./= (sum(fl_heading2.^2).^0.5)
    nir_heading2 = nir_heading .- mean(nir_heading)
    nir_heading2 ./= (sum(nir_heading2.^2).^0.5);
    
    test_xcor = abs.(xcorr(fl_heading2, nir_heading2)) # abs added to correct for 180 flip 
    max_fl_nir_cor, max_fl_nir_cor_idx = findmax(test_xcor)
    
    synchrony_offset = max_fl_nir_cor_idx - length(nir_heading2);
    synchrony_offset_fl = synchrony_offset * 2
    synchrony_offset_sweep = round(Int64, synchrony_offset_fl / 50)
    return synchrony_offset_sweep, length(heading)
end

function keep_valid(x)
    return x[.!isnan.(x)]
end

py"""
from sklearn.linear_model import Ridge

def deconvolve(activity, kernel,alpha=0.5):
    len_kernel = len(kernel)
    kernel_matrix = np.zeros((len(activity), len(activity)+len_kernel-1))
    for i in range(kernel_matrix.shape[0]):
        kernel_matrix[i,i:(i+len_kernel)] = kernel
        
        
    X = kernel_matrix[~np.isnan(activity)]
    y = activity[~np.isnan(activity)]
    
    
    contribution = (kernel_matrix!=0)
    contribution[np.isnan(activity)] = 0
    contribution = np.sum(contribution, axis=0)
    
    clf = Ridge(alpha=alpha,fit_intercept=False,solver="lbfgs",positive=True).fit(X, y)
    activity_deconvolved = clf.coef_
    
    
    return activity_deconvolved, contribution
"""

py"""
import numpy as np
from numba import njit

@njit
def window_baseline(fluoTrace, min_w, p=10):
    numFrames=len(fluoTrace)
    smoothBaseline=np.zeros_like(fluoTrace)
    for i in range(numFrames):
        idx1 = max(0, i - min_w)
        idx2 = min(numFrames-1, i + min_w)
        smoothBaseline[i]=np.nanpercentile(fluoTrace[idx1:idx2],p)
    return smoothBaseline
"""


function fill_nan(vec, desired_length)
    original_length = length(vec)
    new_vec = fill!(Array{Float64,1}(undef, desired_length), NaN)
    new_vec[1:original_length] = vec
    return new_vec
end