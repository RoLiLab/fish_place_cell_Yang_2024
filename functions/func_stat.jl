function ratio_valid(A)
    return sum(.~isnan.(A))/prod(size(A))
end

# shuffle the data, shift
# divide in chunks, then circularly shift
py"""
import numpy as np

def chunk_shuffle_vector(vec, n):
    each_len = int(np.ceil(len(vec)/n))
    chunk_idx_perm = np.random.permutation(n)
    return np.concatenate([vec[each_len*i:each_len*(i+1)] for i in chunk_idx_perm])

def shiftshuffle_vector(vec, at_least = 100):
    y = len(vec)
    perm_idx = np.roll(np.arange(y),np.random.choice(np.arange(at_least,y-at_least)))
    return vec[perm_idx]
"""

function corr_nan(trace1, trace2; at_least_overlap = 2)
    correlation = NaN
    valid_index = (.!isnan.(trace1)).*(.!isnan.(trace2)) #find the overlapping pixel indices
    if sum(valid_index) >= at_least_overlap
        correlation = cor(trace1[valid_index],trace2[valid_index]) #calculate correlation
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

function corr_nan_matrix_2(matrix1, matrix2; at_least_overlap = 2)
    nr_variable_1 = size(matrix1, 2)
    nr_variable_2 = size(matrix2, 2)
    corr_matrix = zeros(Float32, nr_variable_1,nr_variable_2)
    i_j_combo = [(i,j) for i in 1:nr_variable_1 for j in 1:nr_variable_2]
    for (i,j) in i_j_combo
        vector_i = matrix1[:, i]
        vector_j = matrix2[:, j]
        corr_value = corr_nan(vector_i, vector_j; at_least_overlap = 2)
        corr_matrix[i,j] = corr_value
    end
    return corr_matrix
end

function whether_in(vector, collection)
    return [x in collection for x in vector]
end

function sweep_mean(vector::AbstractVector; window_size=ratio_sampling)
    N = length(vector)
    new_N = floor(Int64, N / window_size)
    vector_sweep = reshape(vector[1:new_N*window_size], (window_size, new_N,))
    vector_sweep_mean = dropdims(nanmean(vector_sweep, dims = 1), dims = 1)
    return vector_sweep_mean
end

function sweep_max(vector::AbstractVector; window_size=ratio_sampling)
    N = length(vector)
    new_N = floor(Int64, N / window_size)
    vector_sweep = reshape(vector[1:new_N*window_size], (window_size, new_N,))
    vector_sweep_mean = dropdims(nanmaximum(vector_sweep, dims = 1), dims = 1)
    return vector_sweep_mean
end

function clean_close(when_visit, interval)

    when_visit_real= [] 
    when_old = 1
    for i in 1:length(when_visit) 
        when_new = when_visit[i]
        if i != 1 
            when_old = when_visit[i-1]
        end
        if (when_new - when_old) > interval 
            append!(when_visit_real, when_new) 
        end
    end
    return convert(Vector{Int32}, when_visit_real)
end

function calculate_collect_activity(when_visit_real, how_long, activity_vec)
    collect_activity = fill(NaN32,length(when_visit_real), 2*how_long+1)
    n_frames = size(activity_vec, 1)
    collect_time = []
    collect_time_original = []
    prev_end_time = 0
    for index in 1:length(when_visit_real)
        when = when_visit_real[index]
        if index<length(when_visit_real)
            next_when = when_visit_real[index+1]
            end_time = minimum([when+how_long, n_frames, next_when-40])
        else
            end_time = minimum([when+how_long, n_frames])
        end

        start_time = maximum([when-how_long, 1, prev_end_time+1])
        length_time = end_time-start_time+1
        collect_activity[index, start_time-when+how_long+1:end_time-when+how_long+1] = (activity_vec[start_time:end_time])
        append!(collect_time, [start_time:end_time])
        append!(collect_time_original, [maximum([when-how_long, 1]):minimum([when+how_long, n_frames])])
        prev_end_time = copy(end_time)
    end
    return collect_activity, collect_time,collect_time_original
end
