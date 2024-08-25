module MAP
using PyCall,Statistics,ImageFiltering, NaNStatistics
@pyimport numpy
@pyimport scipy.stats as stats
@pyimport scipy.spatial as spatial

    function calculate_map_direct(activity, which_loc, n_pos; at_least_visit = 1, use_gaussian_filter=false, sigma=1, filter_mask = nothing)
        """Calculate the place map, at the same time returns the matrix (for each bin) for the number of visits and the martix (for each bin) for the summed atcivity.

        Parameters
        ----------
        activity : 1d vector
            neural activity
        which_loc : 1d vector
            which location (which bin) is the fish in
        n_pos : Int
            the number of bins we use along each axis
        at_least_visit : Int (optional), default = 1
            the least number of visits in each bin
        use_gaussian_filter : bool (optional), default = false
            whether or not to use gaussian filter        
        sigma : Float (optional), default = 1
            whether or not to use gaussian filter 
        filter_mask : booolean 2d matrix: n_pos*n_pos (optional), default = nothing
            the mask we use to constrain the filter. E.g. when there are boudaries in the chamber, we don't put any weight on the regions outsides the boundary. The matrix should
            the invalid regions.
    
        Returns
        -------
        cur_map : 2d martix
            place map
        mask_map : 2d martix
            matrix (for each bin) for the number of visits
        activity_num_map : 2d martix
            martix (for each bin) for the summed atcivity

        """
        activity_num_map = calculate_activity_map_digital(activity, which_loc,n_pos; use_gaussian_filter=use_gaussian_filter, sigma=sigma, filter_mask = filter_mask) #the martix (for each bin) for the summed atcivity
        mask_map =calculate_mask_map_digital(which_loc,n_pos; use_gaussian_filter=use_gaussian_filter, sigma=sigma, filter_mask = filter_mask) #the matrix (for each bin) for the number of visits
        cur_map= activity_num_map./mask_map #place maps
        cur_map[mask_map .< at_least_visit] .= NaN; #bins with number of visits smaller than the "at_least_visit" will be set to NaN 
        return cur_map, mask_map, activity_num_map
    end


    function calculate_activity_map_digital(activity::Vector, which_loc::Vector,n_pos; use_gaussian_filter=false, sigma=1, filter_mask = nothing)
       """Calculate the martix (for each bin) for the summed atcivity.

        Parameters
        ----------
        activity : 1d vector
            neural activity
        which_loc : 1d vector
            which location (which bin) is the fish in
        n_pos : Int
            the number of bins we use along each axis
        at_least_visit : Int (optional), default = 1
            the least number of visits in each bin
        use_gaussian_filter : bool (optional), default = false
            whether or not to use gaussian filter        
        sigma : Float (optional), default = 1
            whether or not to use gaussian filter 
        filter_mask : booolean 2d matrix: n_pos*n_pos (optional), default = nothing
            the mask we use to constrain the filter. E.g. when there are boudaries in the chamber, we don't put any weight on the regions outsides the boundary. The matrix should
            the invalid regions.

        Returns
        -------
        P_weighted : 2d martix
            martix (for each bin) for the summed atcivity

        """
        P_weighted = zeros(Float32, n_pos*n_pos) #initialize the output
        bin_count_weighted = numpy.bincount(which_loc, weights=activity) #count the number of visits in each bin
        P_weighted[1:(length(bin_count_weighted)-1)] = bin_count_weighted[2:end];#numpy.bincount starts from 0, but which_loc (indicating which bin) starts from 1
        P_weighted = reshape(P_weighted, n_pos, n_pos); #reshape the output to 2d
        if !isnothing(filter_mask)
            P_weighted[filter_mask].=NaN #invalid regions are set to NaN
        end
        if use_gaussian_filter
            P_weighted_filtered = imfilter(P_weighted, Kernel.gaussian(sigma), NA()) #NA constrained filtering, invalid regions contributes 0 weight
            if !isnothing(filter_mask)
                P_weighted_filtered[filter_mask].=NaN
            end
            P_weighted_filtered = nansum(P_weighted).*P_weighted_filtered./nansum(P_weighted_filtered)
           return P_weighted_filtered
        else
            return P_weighted
        end
    end

    function calculate_mask_map_digital(which_loc::Vector,n_pos; use_gaussian_filter=false, sigma=1, filter_mask = nothing)
       """Calculate the matrix (for each bin) for the number of visits, see "calculate_activity_map_digital" for detailed comments within the function.

        Parameters
        ----------
        which_loc : 1d vector
            which location (which bin) is the fish in
        n_pos : Int
            the number of bins we use along each axis
        at_least_visit : Int (optional), default = 1
            the least number of visits in each bin
        use_gaussian_filter : bool (optional), default = false
            whether or not to use gaussian filter        
        sigma : Float (optional), default = 1
            whether or not to use gaussian filter 
        filter_mask : booolean 2d matrix: n_pos*n_pos (optional), default = nothing
            the mask we use to constrain the filter. E.g. when there are boudaries in the chamber, we don't put any weight on the regions outsides the boundary. The matrix should
            the invalid regions.

        Returns
        -------
        P : 2d martix
            matrix (for each bin) for the number of visits

        """
        P = zeros(Float32, n_pos*n_pos)
        bin_count = Float32.(numpy.bincount(which_loc))    
        P[1:(length(bin_count)-1)] = bin_count[2:end];    
        P = reshape(P, n_pos, n_pos);
        if !isnothing(filter_mask)
            P[filter_mask].=NaN
        end
        if use_gaussian_filter
            P_filtered = imfilter(P, Kernel.gaussian(sigma), NA())
            if !isnothing(filter_mask)
                P_filtered[filter_mask].=NaN
            end
            P_filtered = nansum(P).*P_filtered./nansum(P_filtered)
            return P_filtered
        else
            return P
        end
        return P
    end

    function valid_activity_loc(neural_activity_original::Vector{Float32}, bool_index, loc_digital::Vector{Int64})
       """Find the valid neural acivity and corresponding

        Parameters
        ----------
        neural_activity_original : 1d vector
            original neural acivity (for a single neuron)
        bool_index : BitVector
            specify the range of frames you want to calculate
        loc_digital : 1d vector
            original location (bin) the fish is in    

        Returns
        -------
        P : 2d martix
            matrix (for each bin) for the number of visits

        """
        valid_time_this_neuron = .!(isnan.(neural_activity_original));
        valid_this_neuron = findall(valid_time_this_neuron.*bool_index);
        neural_activity = neural_activity_original[valid_this_neuron]
        which_loc = loc_digital[valid_this_neuron]
        return neural_activity, which_loc
    end

    function map_entropy(map)
        """
        map: 2d matrix, activity map, mean activity in each bin

        """
        # valid_index = findall(.~isnan.(test_map))
        valid_index = findall(map.>=0)
        test_map = map[valid_index]
        entropy_map = NaN
        if length(test_map) !=0
            entropy_map = stats.entropy(test_map, base=length(test_map))
        end
        return entropy_map
    end
    function spatial_info(map, place_mask)
        """
        map: 2d matrix, activity map, mean activity in each bin
        place_mask: 2d matrix, number of visits in each bin

        """
        # valid_index = findall(.~isnan.(test_map))
        valid_index = findall(map.>=0)
        test_map = map[valid_index]
        occupancy = place_mask[valid_index]
        info_content = NaN
        mean_firing = NaN
        if length(test_map) !=0
            visit_ratio = occupancy/sum(occupancy)
            mean_firing = numpy.sum(visit_ratio.*test_map)
            info_content = numpy.nansum(visit_ratio .* test_map .* numpy.log(test_map./mean_firing)./ numpy.log(2), axis=0)
        end
        return info_content, mean_firing
    end
    # function spatial_info(test_map, occupancy)
    #     """
    #     map: 2d matrix, activity map, mean activity in each bin
    #     place_mask: 2d matrix, number of visits in each bin

    #     """
    #     # valid_index = findall(.~isnan.(test_map))
    #     # valid_index = findall(map.>=0)
    #     # test_map = map
    #     # occupancy = place_mask
    #     info_content = NaN
    #     if sum(test_map.>0) > 10
    #         visit_ratio = occupancy/numpy.nansum(occupancy)
    #         mean_firing = numpy.nansum(visit_ratio.*test_map)
    #         info_content = numpy.nansum(visit_ratio .* test_map .* numpy.log(test_map./mean_firing)./ numpy.log(2))
    #     end
    #     return info_content
    # end

    function corr_peak_fast(map)
        len_x = size(map, 1)
        len_y = size(map, 2)
        dx_dy_all = [(0,1), (1,0), (0,-1), (-1,0)]
        #initialize the output, reformat the map
        auto_correlate2d = zeros(Float32, 4)
        map1_new = fill(NaN32, 3*len_x-2, 3*len_y-2)
        map1_new[len_x:2*len_x-1, len_y:2*len_y-1] .=map;
        i = 1
        for (dx, dy) in dx_dy_all
            map2_new = fill(NaN32, 3*len_x-2, 3*len_y-2)
            map2_new[len_x+dx:2*len_x-1+dx, len_y+dy:2*len_y-1+dy] .=map
            valid_index = findall((.!isnan.(map1_new)).*(.!isnan.(map2_new))) #find the overlapping pixel indices
            auto_correlate2d[i] = numpy.corrcoef(map1_new[valid_index],map2_new[valid_index])[1,2] #calculate correlation
            i += 1
        end
        return mean(auto_correlate2d)
    end
end


# function digital2map(which_loc::Vector,n_pos::Int)
#     P = zeros(Float32, n_pos*n_pos)
#     bin_count = Float32.(numpy.bincount(which_loc))    
#     P[1:(length(bin_count)-1)] = bin_count[2:end];
#     return digital2map


function spatial_variance(map)
    """
    map: 2d matrix, activity map, mean activity in each bin
    
    """
    test_map = copy(map)
    if sum(isfinite.(test_map)) == 0
        return NaN32
    end
    test_map .-= minimum(test_map[isfinite.(test_map)]);

    test_map_xy = findall(isfinite.(test_map)); 
    

    test_map_x = [test_map_xy[x][1] for x = 1:length(test_map_xy)];
    test_map_y = [test_map_xy[x][2] for x = 1:length(test_map_xy)];

    test_map_p = test_map ./ nansum(test_map); 
    test_map_p = test_map_p[test_map_xy];


    test_map_mu_x = sum(test_map_x .* test_map_p)
    test_map_mu_y = sum(test_map_y .* test_map_p)

    map_var = sum(((test_map_x .- test_map_mu_x).^2 + (test_map_y .- test_map_mu_y).^2) .* test_map_p)
    return map_var
end


# need to change to chunks only on valid indices
function chunk_maps(neural_activity_original, loc_digital, bool_index_all; at_least_visit = 1, use_gaussian_filter=false, sigma=1, filter_mask=nothing)
    how_many_chunk = length(bool_index_all)
    map_all = fill(NaN32, n_pos, n_pos, how_many_chunk)
    for i in 1:how_many_chunk
        bool_index = bool_index_all[i]
        neural_activity, which_loc = MAP.valid_activity_loc(neural_activity_original, bool_index,loc_digital)
        map_all[: ,:, i],_,_  = MAP.calculate_map_direct(neural_activity, which_loc,n_pos; at_least_visit = at_least_visit, use_gaussian_filter=use_gaussian_filter, sigma=sigma, filter_mask=filter_mask)
    end
    return map_all
end

function generate_chunks(how_many_chunk=2; n_sweeps=n_sweeps,filter=nothing)
    each_chunk_length = floor(Int64,n_sweeps/how_many_chunk)
    bool_index_all = []
    for i in 1:how_many_chunk
        bool_index_one = falses(n_sweeps)
        bool_index_one[1+each_chunk_length*(i-1):each_chunk_length*i].=true
        if .!isnothing(filter)
            bool_index_one[filter].=false
        end
        append!(bool_index_all, [bool_index_one])
    end
    return bool_index_all
end

function generate_chunks_fine(how_many_chunk=2; n_sweeps=n_sweeps, filter=nothing)
    each_chunk_length = floor(Int64,n_sweeps/how_many_chunk)
    chunk_label = zeros(Int32, n_sweeps)
    for i in 1:how_many_chunk
        chunk_label[1+each_chunk_length*(i-1):each_chunk_length*i].=i
    end
    set_1_label = sort(numpy.random.choice(1:how_many_chunk,Int(how_many_chunk/2),replace=false))
    bool_index_1 = falses(n_sweeps)
    bool_index_2 = falses(n_sweeps)
    for i in 1:how_many_chunk
        if i in set_1_label
            bool_index_1[1+each_chunk_length*(i-1):each_chunk_length*i].= true
        else
            bool_index_2[1+each_chunk_length*(i-1):each_chunk_length*i].= true
        end
    end
    if .!isnothing(filter)
        bool_index_1[filter].=false
        bool_index_2[filter].=false
    end
    return [bool_index_1, bool_index_2]
end

function corr_2d(map1, map2; threshold =10)
    auto_correlate2d = NaN32
    valid_index = findall((.!isnan.(map1)).*(.!isnan.(map2))) #find the overlapping pixel indices
    if length(valid_index) >= threshold
        auto_correlate2d = numpy.corrcoef(map1[valid_index],map2[valid_index])[1,2] #calculate correlation
    end
    return auto_correlate2d
end

function corr_2d_original(map1, map2; threshold =10, mask = nothing)
    corr_maps = NaN32
    valid_index = findall((.!isnan.(map1)).*(.!isnan.(map2))) #find the overlapping pixel indices
    if !isnothing(mask)
        valid_index = findall((.!isnan.(map1)).*(.!isnan.(map2)).*(.!isnan.(mask)))
    end
    if length(valid_index) >= threshold
        cov_maps = mean((map1[valid_index] .- numpy.nanmean(map1)).*(map2[valid_index] .- numpy.nanmean(map2)))
        var_maps = sqrt(numpy.nanmean((map1 .- numpy.nanmean(map1)).^2))*sqrt(numpy.nanmean((map2 .- numpy.nanmean(map2)).^2))
        corr_maps = cov_maps/var_maps 
    end
    return corr_maps
end

@pyimport numpy
function autocorr2d(map1, map2, threshold =30)
    """
    map1 & map2: 2d matrix, activity map, mean activity in each bin
    
    """
    #the ranges we use to shift the map
    len_x = size(map1, 1)
    len_y = size(map1, 2)
    xrange = -(len_x-1):(len_x-1)
    yrange = -(len_y-1):(len_y-1)
    #initialize the output, reformat the map
    auto_correlate2d = zeros(Float64, length(xrange), length(yrange))
    map1_new = fill(NaN32, 3*len_x-2, 3*len_y-2)
    map1_new[len_x:2*len_x-1, len_y:2*len_y-1] .=map1;
    for (i, dx) in enumerate(xrange)
        for (j, dy) in enumerate(yrange)
            map2_new = fill(NaN32, 3*len_x-2, 3*len_y-2)
            map2_new[len_x+dx:2*len_x-1+dx, len_y+dy:2*len_y-1+dy] .=map2
            valid_index = findall((.!isnan.(map1_new)).*(.!isnan.(map2_new))) #find the overlapping pixel indices
            if length(valid_index)>threshold
                auto_correlate2d[i,j] = numpy.corrcoef(map1_new[valid_index],map2_new[valid_index])[1,2] #calculate correlation
            end
        end
    end
    return auto_correlate2d
end

# @pyimport scipy.signal as signal
function corr_peak_half_decay(auto_correlate2d, dr_reshaped)
    
    auto_correlate2d_reshaped = reshape(auto_correlate2d, length(auto_correlate2d));
    p = sortperm(dr_reshaped)
    dr_sorted = dr_reshaped[p]
    auto_correlate2d_sorted = auto_correlate2d_reshaped[p]
    dr_unique = numpy.unique(dr_sorted);

    auto_correlate2d_averaged = [mean(auto_correlate2d_sorted[dr_sorted.==dr_one]) for dr_one in dr_unique];
    # figure()
    # plot(dr_unique,auto_correlate2d_averaged )
    max_corr = maximum(auto_correlate2d_averaged[2:end])
    # max_corr = auto_correlate2d_averaged[2]
    
    candicate_r = collect(0:0.01:maximum(dr_unique))
    half_max = max_corr/2
    sem = [sum((auto_correlate2d_averaged.-half_max).*((dr_unique.>r_c).*2 .-1)) for r_c in candicate_r];
    sem_minimum_index = findall(sem.==minimum(sem))
    half_decay_distance = mean(candicate_r[sem_minimum_index])
    
    return max_corr, half_decay_distance
end

function map_stability(maps)
    nr_map = size(maps,3)
    corr_all = fill(NaN32, Int(nr_map*(nr_map-1)/2))
    t = 1
    for i in 1:(nr_map-1)
        for j in i+1:nr_map
            corr_all[t] = corr_2d(maps[:,:,i], maps[:,:,j])
            t+=1
        end
    end
    return nanmean(corr_all)
end
@pyimport networkx as nx
function connected_cluster_label(matrix)
    if typeof(matrix) in [Array{Float32, 2}, Array{Float64, 2}]
        nr_node = size(matrix, 1)
        G = nx.from_numpy_matrix(matrix)
    else
        nr_node = matrix.shape[1]
        G = nx.from_scipy_sparse_array(matrix)
    end
    # for i in 1:nr_node
    #     matrix[i,i] = false
    # end
    # 
    
    G_components = [py"list"(c).+1 for c in nx.connected_components(G)]
    component_label = zeros(Int32, nr_node)
    # random_labels = numpy.random.permutation(1:length(G_components))
    random_labels = 1:length(G_components)
    for i in 1:length(G_components)
        component_label[G_components[i]].=random_labels[i]
    end
    return component_label
end

function corr_multiple(array)
    # the same as nancor?
    # nr_feature * nr_sample
    nr_feature = size(array, 1)
    nr_sample  = size(array, 2)
    matrix = fill(NaN32, nr_sample, nr_sample)
    @showprogress for i in 1:nr_sample
        for j in 1:nr_sample
            matrix[i, j] = corr_2d(array[:, i], array[:, j])
        end
    end
    return matrix
end


# function find_same_neuron(cell_index, X_all, Y_all, Z_all, A_dF)
#     place_cell_xy = hcat(X_all[cell_index], Y_all[cell_index])
#     place_cell_z = hcat(Z_all[cell_index])
#     # corr_activity = corr_multiple(A_dF[:,cell_index])
#     corr_activity = corr_nan_matrix(A_dF[:,cell_index])
#     same_loc_matrix = (spatial.distance_matrix(place_cell_xy,place_cell_xy) .<2) .* (spatial.distance_matrix(place_cell_z,place_cell_z) .<2)
#     same_activity_matrix = (corr_activity.>0.7);
#     same_neuron_matrix = same_loc_matrix.*same_activity_matrix
#     component_label = connected_cluster_label(same_neuron_matrix);
#     return component_label
# end

@pyimport sklearn.neighbors as neighbors
@pyimport scipy.sparse as sparse
function find_same_neuron(loc_x, loc_y, loc_z, activity;corr_thres = 0.7)
    place_cell_xy = hcat(loc_x, loc_y)
    place_cell_z = hcat(loc_z)
    #calculate spatial closeness first
    xy_matrix = neighbors.radius_neighbors_graph(place_cell_xy, 1.42, mode="connectivity",
                               include_self=false).astype(numpy.bool_)
    z_matrix = neighbors.radius_neighbors_graph(place_cell_z, 1.42, mode="connectivity",
                               include_self=false).astype(numpy.bool_)

    same_loc_matrix = xy_matrix.multiply(z_matrix)
    
    sparse_x,sparse_y,_ = sparse.find(sparse.triu(same_loc_matrix, k=1));
    sparse_x = sparse_x.+1
    sparse_y = sparse_y.+1 #correct for the python indexing;
    # now only compute correlation for those who are close
    same_neuron_matrix = same_loc_matrix.copy()
    for i in 1:length(sparse_x)
        x = sparse_x[i]
        y = sparse_y[i]
        corr_value = corr_nan(activity[:,x], activity[:,y]; at_least_overlap = 1)
        # print(x, y, corr_value)
        if !(corr_value.>corr_thres) 
            same_neuron_matrix[x, y] = false
            same_neuron_matrix[y, x] = false
        end
    end
    same_neuron_matrix.eliminate_zeros()
    component_label = connected_cluster_label(same_neuron_matrix);
    return component_label
end

function divide_long(how_long, max_length)
    nr_piece = ceil(Int32, how_long/max_length)
    return [((i-1)*max_length+1):minimum([((i)*max_length), how_long]) for i in 1:nr_piece]
end

@pyimport itertools

function divide_long_corr(which_roi, node_z; max_length=5,A_dF=A_dF)
    nr_cuts = ceil(Int32,length(which_roi)/max_length)-1
    p = sortperm(node_z)
    corr_matrix = corr_nan_matrix(A_dF[:,which_roi[p]])
    all_possible_cuts = []
    where_cut_all = []
    for combo in itertools.combinations(1:(length(which_roi)-1), nr_cuts)
        where_cut = reduce(vcat, [0, [combo[i] for i in 1:length(combo)],length(which_roi)])
        cell_lengths = diff(where_cut)
        if maximum(cell_lengths)<=max_length
            append!(all_possible_cuts, [combo])
            append!(where_cut_all, [where_cut])
        end
    end
    clustering_index_all = []
    for where_cut in where_cut_all
        clustering_matrix_within = fill(NaN32, size(corr_matrix));
        clustering_matrix_cross = copy(corr_matrix);
        clustering_index = 0
        for i in 1:(length(where_cut)-1)
            clustering_matrix_cross[where_cut[i]+1:where_cut[i+1],where_cut[i]+1:where_cut[i+1]].=NaN
            clustering_matrix_within[where_cut[i]+1:where_cut[i+1],where_cut[i]+1:where_cut[i+1]].=corr_matrix[where_cut[i]+1:where_cut[i+1],where_cut[i]+1:where_cut[i+1]]
            # clustering_index+=mean(clustering_matrix[where_cut[i]+1:where_cut[i+1],where_cut[i]+1:where_cut[i+1]])
            # clustering_matrix[where_cut[i]+1:where_cut[i+1],where_cut[i]+1:where_cut[i+1]] .=  -clustering_matrix[where_cut[i]+1:where_cut[i+1],where_cut[i]+1:where_cut[i+1]]
        end
        clustering_index = numpy.nanmean(clustering_matrix_within)-numpy.nanmean(clustering_matrix_cross)
        append!(clustering_index_all, clustering_index)
    end
    best_cut = all_possible_cuts[findmax(clustering_index_all)[2]]
    best_where_cut = where_cut_all[findmax(clustering_index_all)[2]]
    cell_lengths = diff(where_cut_all[findmax(clustering_index_all)[2]]);
    cell_index = [best_where_cut[i]+1:best_where_cut[i+1] for i in 1:(length(best_where_cut)-1)]
    
    return cell_index
end

function deal_long_neuron(component_label, z_loc;valid_neurons=valid_neurons, max_length = 5,A_dF=A_dF)
    
    labels, counts = numpy.unique(component_label, return_counts = true);

    new_label = 1
    which_too_large = labels[findall(counts .> max_length)]
    component_label_small = fill(NaN32, length(component_label))
    for (i, c_label) in enumerate(labels)
        c_size = counts[c_label]
        which_nodes = findall(component_label.==c_label)
        if c_size > max_length
            node_z = z_loc[which_nodes]
            p = sortperm(node_z)
            # indices_ = divide_long(c_size, max_length)
            indices_ = divide_long_corr(valid_neurons[which_nodes], node_z;A_dF=A_dF)
            for index in indices_
                # println(new_label)
                component_label_small[which_nodes[p[index]]] .= new_label
                new_label +=1
            end
        else
            component_label_small[which_nodes] .= new_label
            new_label +=1
        end
    end
    component_label_small = Int32.(component_label_small);
    
    random_labels = 1:maximum(component_label_small)
    # random_labels = numpy.random.permutation(1:maximum(component_label_small))
    component_label_small_random = fill(NaN32, length(component_label_small))
    for i in 1:maximum(component_label_small)
        component_label_small_random[component_label_small.==i].=random_labels[i]
    end
    component_label_small_random = Int32.(component_label_small_random);
    return component_label_small_random
end

function normalizing_map(map, min_zero = true)
    if min_zero
        min_map = numpy.nanmin(map)
    else
        min_map = 0
    end
    normalized_map = map .- min_map
    sum_map = nansum(normalized_map)
    if sum_map != 0
        normalized_map = normalized_map./sum_map
        return normalized_map
    else 
        return []
    end
end


function space_by_weight(weights, min_bin, max_bin)
    weights_norm = weights./sum(weights)
    whole_range = max_bin-min_bin
    nr_bins = length(weights)
    bins=[min_bin]
    bins = Float32.(bins)

    for i in 1:(nr_bins-1)
        weight= weights_norm[i]
        append!(bins, bins[i]+whole_range*weight)
    end
    bins = vcat(bins, max_bin)
    return bins
end

function ego_peak_map(which_neuron, which_indices, ax; nr_r_bins=10,  nr_a_bins=24, at_least_visit=5, r_spacing="weight")
    map = neuron_map_all_original[:,:,which_neuron]
    # imshow(map[:,1:8], origin="lower")
    peak_loc_idx = findall(map.==numpy.nanmax(map))[1]
    peak_loc = [x_bins_mid[peak_loc_idx[1]], y_bins_mid[peak_loc_idx[2]]]

    fish_distance_peak = distance_from(hcat(x_fish_sweep_mean, y_fish_sweep_mean), peak_loc)
    fish_angle_peak = angle_from(hcat(x_fish_sweep_mean, y_fish_sweep_mean), peak_loc) .+pi;

    facing = heading_sweep_mean.-fish_angle_peak
    facing[facing.>pi] .-= 2*pi
    facing[facing.<-pi] .+= 2*pi

    valid_indices = intersect(findall(.!isnan.(A_dF[:,which_neuron])), 3601:n_sweeps)
    azimut = facing[valid_indices]
    azimut[azimut.<0].+=2*pi
    radius = fish_distance_peak[valid_indices]
    # define binning
    if r_spacing == "linear"
        rbins = numpy.linspace(0,maximum(radius), nr_r_bins)
    else
        rbins = space_by_weight(1 ./numpy.linspace(1,nr_r_bins, nr_r_bins), 0,maximum(radius))
    end
    abins = numpy.linspace(0,2*pi, nr_a_bins)

    #calculate histogram
    hist_count, _, _ = numpy.histogram2d(azimut, radius, [abins, rbins])
    hist_count[hist_count.<at_least_visit].=NaN
    hist_activity, _, _ = numpy.histogram2d(azimut, radius, [abins, rbins], weights = A_dF[valid_indices,which_neuron])
    hist_ = hist_activity./hist_count
    A, R = numpy.meshgrid(abins, rbins)

    # plot
    pc = ax.pcolormesh(A, R, hist_', cmap="viridis")
    # clb = colorbar(pc)
    # clb.ax.set_title("dF")
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


function find_peak(map)
    """
    map: 2d matrix, activity map, mean activity in each bin
    
    """
    test_map = copy(map)
    map_max = nanmaximum(test_map)
    peak_loc = findall(test_map.==map_max)[1]
    return peak_loc[1], peak_loc[2]
end

"""
Functions for detecting the edges
"""
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


function map_components(map; coarse=true, coarse_sigma=1, threshold=95)
    """
    find out how many components are in the map, thresholded by 95 percentile
    
    Parameters
    ----------
    coarse : whether or not to use gaussian filter
    coarse_sigma : gaussian filter size
    
    Returns
    -------
    components_size : size of each component
    img_label : each component is labeled with one number
    coarse_map : gaussian filtered map, if coarse = false, return original image
    """
    coarse_map = copy(map)
    if coarse
        coarse_map = imfilter(coarse_map, Kernel.gaussian(coarse_sigma), NA())
    end
    coarse_map[isnan.(map)].=NaN
    map_max = numpy.nanmax(coarse_map)
    thresholded_map = copy(coarse_map)
    invalid_index = thresholded_map.<nanpctile(reshape(thresholded_map, length(thresholded_map)), threshold)
    coarse_map[invalid_index].=NaN
    thresholded_map[invalid_index].=NaN
    thresholded_map[.!isnan.(thresholded_map)].=1
    thresholded_map[isnan.(thresholded_map)].=0
    img_label = label_components(Int64.(thresholded_map));
    img_bg_label_area = component_lengths(img_label)
    components_size = img_bg_label_area[2:end]
    return components_size, img_label, coarse_map
end


function map_classification(maps; n_clusters=8)
    """
    classify maps based on the peak location, return class and peak locations
    """
    peak_loc_map =  fill(NaN32,2, size(maps,3))
    for i in 1:size(maps,3)
        peak_loc_map[:,i] .= find_com(maps[:,:,i])
    end
    KMeans_features = cluster.KMeans(n_clusters = n_clusters).fit(hcat(peak_loc_map'))
    return KMeans_features.labels_, peak_loc_map
end
    
function find_com(which_map)
    """
    com of a map, input is normally thresholded map
    
    """
    test_map = copy(which_map)
    test_map .-=numpy.nanmin(test_map)

    test_map_xy = findall(isfinite.(test_map)); 

    test_map_x = [test_map_xy[x][1] for x = 1:length(test_map_xy)];
    test_map_y = [test_map_xy[x][2] for x = 1:length(test_map_xy)];

    test_map_p = test_map ./ nansum(test_map); 
    test_map_p = test_map_p[test_map_xy];


    test_map_mu_x = sum(test_map_x .* test_map_p)
    test_map_mu_y = sum(test_map_y .* test_map_p)

    
    return (test_map_mu_x, test_map_mu_y)
end


function count_mode(map; filter=true,coarse_sigma=2, mask=nothing,min_component_size = 5)
    test_map = largest_component(map)
    if filter
        test_map = filter_map(test_map;coarse_sigma=coarse_sigma, mask=mask)
    end
    nr_component = find_components(test_map;min_component_size=min_component_size)
    return nr_component
end

function largest_component(map)
    test_map = copy(map)
    test_map[.!isnan.(test_map)].=1
    test_map[isnan.(test_map)].=0
    img_label = label_components(Int64.(test_map));
    img_bg_label_area = component_lengths(img_label)

    max_compoennt_label = findmax(img_bg_label_area[2:end])[2]
    map_largest_component = copy(map)
    map_largest_component[img_label.!=max_compoennt_label].=NaN
    return map_largest_component
end

function filter_map(map;coarse_sigma=2, mask=nothing)
    map_masked = copy(map)
    if ! isnothing(mask)
        map_masked[.!mask].=NaN
    end
    
    coarse_map = imfilter(map, Kernel.gaussian(coarse_sigma), NA())
    # coarse_map[isnan.(map)].=NaN
    return coarse_map
end
function find_most_common(array)
    values, counts = numpy.unique([1,1,2,2,2,3], return_counts=true)
    ind = numpy.argmax(counts)+1
    return values[ind]
end

function find_components(which_map;min_component_size=5)
    # maps should be gaussian filtered with sigma 1
    components = fill(NaN, 100)
    components_figure = fill(NaN, 100, size(which_map,1), size(which_map,2))
    for p in 50:100
        test_map = copy(which_map)
        test_map[test_map.<nanpctile(reshape(test_map, length(test_map)), p)].=NaN
        test_map[.!isnan.(test_map)].=1
        test_map[isnan.(test_map)].=0
        img_label = label_components(Int64.(test_map));

        img_bg_label_area = component_lengths(img_label)

        real_components = (1:maximum(img_label))[img_bg_label_area[2:end].>=min_component_size]
        # components_value = [mean(which_map[img_label.==which_component]) for which_component in real_components]
        components_value = [numpy.nanmean(which_map[img_label.==which_component]) for which_component in real_components]

        if length(components_value)>1
            components[p] = sum((components_value.>maximum(components_value) *4/5))
            components_figure[p,:,:] = keep_index(img_label, real_components[components_value.>maximum(components_value) *4/5])

        elseif length(components_value)==1
            components[p] =1
            components_figure[p,:,:] = keep_index(img_label, real_components)
        end
    end
    
    components_median = round(Int32, numpy.nanmedian(components[components.>0]))
    
    sum_components = numpy.sum(components_figure[findall(components.==components_median), :,:], axis=0)
    sum_components[.!isnan.(sum_components)].=1
    sum_components[isnan.(sum_components)].=0
    img_label = label_components(Int64.(sum_components));
    img_bg_label_area = component_lengths(img_label)
    real_components = (1:maximum(img_label))
    components_value = [mean(which_map[img_label.==which_component]) for which_component in real_components]
    components_loc = [component_loc(img_label, which_component) for which_component in real_components]
    return components_median, components_value, components_loc
end

function component_loc(img_label, which_component) 
    component_index = findall(img_label.==which_component)
    mean_x = mean([component_index[i][1] for i in 1:length(component_index)])
    mean_y = mean([component_index[i][2] for i in 1:length(component_index)])
    return (mean_x, mean_y)
end

function keep_index(img_label, which_components)
    components_list = 1:maximum(img_label)
    img_label_new = Float32.(copy(img_label))
    for i in components_list
        if !(i in which_components)
            img_label_new[img_label_new.==i].=0
        end
    end
    img_label_new[img_label_new.==0].=NaN
    # img_label_new[img_label_new.!=0].=1
    return img_label_new
end

@pyimport scipy.ndimage as ndimage

function solve_geometry_old(x_fish, y_fish, C, img_bg;method="trajectory",background_threshold_low=0,background_threshold_high=250)
    if method == "trajectory"
        xy_projection = zeros(size(img_bg[:,:,end])); 

        for i = findall(C.>0.4)
        # for i = 1:length(x_fish)

            x = round.(Int64, x_fish[i])
            y = round.(Int64, y_fish[i])
            xy_projection[x, y] = 1;
        end

        σ = 3.0f0
        kern_σ = KernelFactors.gaussian((σ, σ))
        img_xy = imfilter(xy_projection, kern_σ);
        img_xy = (img_xy.>0);
    else
        img_bg_end = img_bg[:,:,end];
        l = size(img_bg_end)[1]
        w = size(img_bg_end)[2]


        figure()
        img_hist(img_bg_end, th_l=0, th_h=500, bins=200);


        #assumption 1: chamber size larger than 100000
        #assumption 2: no attachment to the edge
        #assumption 3: center closest the image center


        minimim_component_size = 10000

        img_xy = (img_bg_end .< background_threshold_high) .* (img_bg_end .> background_threshold_low); 
        img_xy .= erode(img_xy);
    end
    figure()
    imshow(img_xy)
    img_label = label_components(img_xy);
    img_bg_label_area = component_lengths(img_label)
    candidate_roi_idx = (findall(img_bg_label_area .> 100000).-1); #julia1.7
    # candidate_roi_idx = (findall(img_bg_label_area .> 100000));
    edge_covering = zeros(length(candidate_roi_idx))
    for (index, i) in enumerate(candidate_roi_idx)
        roi = (img_label .== i)
        edge_covering[index] = sum(roi[1,:])+sum(roi[end,:])+sum(roi[:,1])+sum(roi[:,end])
    end
    roi_idx = candidate_roi_idx[edge_covering.==minimum(edge_covering)];
    
    println(roi_idx)

    chamber_roi = (img_label .== roi_idx);

    for i = 1:3
        chamber_roi .= dilate(closing(chamber_roi))
        # fig, ax = PyPlot.subplots(1,1)
        # ax.imshow(chamber_roi)
    end
    

    #find countour and center
    countour = find_contours(chamber_roi)[1];
    countour_array = hcat(getindex.(countour, 1), getindex.(countour,2));
    center_loc = mean(countour_array, dims=1);
    BoolBorder = zeros(Int64, size(chamber_roi));
    for ind in countour
        BoolBorder[ind] = 1
    end

    chamber_roi = convert(Array{Int64}, chamber_roi);
    
    chamber_roi = ndimage.binary_fill_holes(BoolBorder)
    chamber_roi = convert(Array{Int64}, chamber_roi)
    
    figure()
    imshow(chamber_roi')
    scatter(center_loc[1], center_loc[2], color="r", s=10)
    
    chamber_roi_xy = findall(chamber_roi.!=0)
    chamber_roi_x = [xy[1] for xy in chamber_roi_xy][1:1000:end]
    chamber_roi_y = [xy[2] for xy in chamber_roi_xy][1:1000:end];
    return countour_array, center_loc,chamber_roi, chamber_roi_x, chamber_roi_y
end

function solve_geometry(x_fish, y_fish, C, img_bg;method="trajectory",background_threshold_low=0,background_threshold_high=250)
    if method == "trajectory"
        xy_projection = zeros(size(img_bg[:,:,end])); 

        for i = findall(C.>0.4)
        # for i = 1:length(x_fish)

            x = round.(Int64, x_fish[i])
            y = round.(Int64, y_fish[i])
            xy_projection[x, y] = 1;
        end

        σ = 3.0f0
        kern_σ = KernelFactors.gaussian((σ, σ))
        img_xy = imfilter(xy_projection, kern_σ);
        img_xy = (img_xy.>0);
    else
        img_bg_end = img_bg[:,:,end];
        l = size(img_bg_end)[1]
        w = size(img_bg_end)[2]


        figure()
        img_hist(img_bg_end, th_l=0, th_h=500, bins=200);


        #assumption 1: chamber size larger than 100000
        #assumption 2: no attachment to the edge
        #assumption 3: center closest the image center


        minimim_component_size = 10000

        img_xy = (img_bg_end .< background_threshold_high) .* (img_bg_end .> background_threshold_low); 
        img_xy .= erode(img_xy);
    end
    figure()
    imshow(img_xy)
    img_label = label_components(img_xy);
    img_bg_label_area = component_lengths(img_label)
    # candidate_roi_idx = (findall(img_bg_label_area .> 100000).-1); #julia1.7
    candidate_roi_idx = (findall(img_bg_label_area .> 100000));
    edge_covering = zeros(length(candidate_roi_idx))
    for (index, i) in enumerate(candidate_roi_idx)
        roi = (img_label .== i)
        edge_covering[index] = sum(roi[1,:])+sum(roi[end,:])+sum(roi[:,1])+sum(roi[:,end])
    end
    roi_idx = candidate_roi_idx[edge_covering.==minimum(edge_covering)];
    
    println(roi_idx)

    chamber_roi = (img_label .== roi_idx);

    for i = 1:3
        chamber_roi .= dilate(closing(chamber_roi))
        # fig, ax = PyPlot.subplots(1,1)
        # ax.imshow(chamber_roi)
    end
    

    #find countour and center
    countour = find_contours(chamber_roi)[1];
    countour_array = hcat(getindex.(countour, 1), getindex.(countour,2));
    center_loc = mean(countour_array, dims=1);
    BoolBorder = zeros(Int64, size(chamber_roi));
    for ind in countour
        BoolBorder[ind] = 1
    end

    chamber_roi = convert(Array{Int64}, chamber_roi);
    
    chamber_roi = ndimage.binary_fill_holes(BoolBorder)
    chamber_roi = convert(Array{Int64}, chamber_roi)
    
    figure()
    imshow(chamber_roi')
    scatter(center_loc[1], center_loc[2], color="r", s=10)
    
    chamber_roi_xy = findall(chamber_roi.!=0)
    chamber_roi_x = [xy[1] for xy in chamber_roi_xy][1:1000:end]
    chamber_roi_y = [xy[2] for xy in chamber_roi_xy][1:1000:end];
    return countour_array, center_loc,chamber_roi, chamber_roi_x, chamber_roi_y
end


# define firing fields
function map_field(which_map; threshold = 8/10, bottom_activity= 0)
    top_activity = numpy.nanpercentile(which_map, 95)
    field_threshold = (top_activity-bottom_activity)*threshold +bottom_activity
    thresholded_map = copy(which_map)
    valid_index = thresholded_map.>(field_threshold)

    return valid_index
end


function dynamic_range(which_map)
    bottom_activity = numpy.nanpercentile(which_map, 5)
    top_activity = numpy.nanpercentile(which_map, 95)

    return (top_activity - bottom_activity)/bottom_activity
end

function map_field_size(which_map; threshold = 8/10)
    valid_index = map_field(which_map; threshold= threshold)
    return sum(valid_index)
    
end

function map_valid_field_size(which_map)
    components_peaks, img_label_valid, valid_components = map_components_peak(which_map)
    return sum(img_label_valid.!=0)
    
end

# define mode peaks
function component_peak(img_label, which_component, which_map) 
    not_component_index = findall(img_label.!=which_component)
    
    map_this_component = copy(which_map)
    map_this_component[not_component_index] .= NaN32

    return find_com(map_this_component)
end


function map_components_peak(which_map; threshold = 8/10, components_size_threshold = 20)
    firing_field = map_field(which_map; threshold = threshold)
    
    
    img_label = label_components(Int64.(firing_field));
    img_bg_label_area = component_lengths(img_label)
    # components_size = img_bg_label_area[2:end]

    components_size = img_bg_label_area[1:end] # fix due to updates

    valid_components = findall(components_size.>components_size_threshold)
    nr_component = length(valid_components)

    img_label_valid = copy(img_label)
    img_label_valid[.! whether_in(img_label, valid_components)] .=0
    
    components_peaks = [component_peak(img_label, which_component, which_map) for which_component in valid_components]
    return components_peaks, img_label_valid, valid_components
end

function map_peak(which_map; threshold = 8/10)
    firing_field = map_field(which_map; threshold = threshold)
    
    map_firing_field = fill(NaN32, size(which_map))
    map_firing_field[firing_field] .= which_map[firing_field]

    return find_com(map_firing_field)
end

function Moran_I(which_map, w, map_pixels)
    map_values = which_map[map_pixels]
    map_values_mean = mean(map_values)
    map_values_var = sum((map_values .- map_values_mean).^2)
    N = length(map_pixels)
    W = sum(w)
    
    I = 0
    for i in 1:length(map_pixels)
        for j in 1:length(map_pixels)
            I += (map_values[i] - map_values_mean)*(map_values[j] - map_values_mean)*w[i,j]
        end
    end

    I = I*N/W/map_values_var
    
    
    return I
end

# to find neurons with confined place fields (return index of place_cell_index) and the best peak of all place_cell_index
function neuron_with_valid_peak(place_cell_index, place_map_all;threshold=8/10)
    best_peak_all = fill(NaN32, 2, length(place_cell_index))
    firing_field_size_ratio_all = fill(NaN32, length(place_cell_index))
    nr_modes = fill(NaN32, length(place_cell_index))
    for (i, which_neuron) in enumerate(place_cell_index)
        mask_valid = isfinite.(place_map_all[:,:,which_neuron])
        which_map = place_map_all[:,:,which_neuron]
        components_peaks, img_label_valid, valid_components = map_components_peak(which_map; threshold = threshold, components_size_threshold = 20)
        firing_field_size_ratio = sum(img_label_valid.!=0)/sum(mask_valid)
        peak_values = [numpy.nanpercentile(which_map[img_label_valid .== valid_components[i]], 95) for i in 1:length(valid_components)]
        if length(peak_values) .> 0
            best_peak = components_peaks[findmax(peak_values)[2]]
            firing_field_size_ratio_all[i] = firing_field_size_ratio
            best_peak_all[:,i] .= best_peak
            nr_modes[i] = length(valid_components)
        end
    end
    # which_neuron = intersect(findall(firing_field_size_ratio_all.<0.3), findall(nr_modes .==1))

    which_neuron = intersect(findall(firing_field_size_ratio_all.<0.3))
    confined_place_cell_index = place_cell_index[which_neuron];
    
    return which_neuron, best_peak_all
end

function whether_in(vector, collection)
    return [x in collection for x in vector]
end
