module Decoder

using PyCall, PyPlot, ProgressMeter, NaNStatistics
using Random, Distributions
using Images
np = pyimport("numpy")
include("../project_place_cell/functions/func_map.jl")


    function rolling_decoder(include_mask, A_dF_place_cells, loc_digital, bins::Int, shift::Int, n_pos::Int; decoding_window=120, space=120)
    #perform the rolling 1min prediction with 1 min buffer on each side 
        
        n_sweeps = size(A_dF_place_cells, 1)

        rolling_predicted = fill(NaN, n_sweeps-shift, 2)
        rolling_valid = fill(NaN, n_sweeps-shift)
        rolling_var = fill(NaN, n_sweeps-shift)


        #shift activity against location
        A_dF_place_cells_shifted, loc_shifted, include_mask_shifted = Decoder.shift_time(shift, A_dF_place_cells, [loc_digital, include_mask])

        for_density_maps = copy(Decoder.make_basis_maps(A_dF_place_cells_shifted, loc_shifted, include_mask_shifted, n_pos, true))
        
        density_map = Float32.(sum([for_density_maps[:,:,i] .> nanpctile(for_density_maps[:,:,i], 95) for i in 1:size(for_density_maps, 3)])) #was 95pctile
        density_map[density_map .== 0] .= 1

        density_map .= density_map .^(1/3)
        
    
        A_dF_place_cells_binned = Decoder.bin_activity(bins, A_dF_place_cells_shifted)
        
        for i in space:decoding_window:n_sweeps-(decoding_window+space)-shift
            train_on = trues(n_sweeps-shift)
            train_on[i-space+1:i+(decoding_window+space)] .= false
            neuron_map = Decoder.make_basis_maps(A_dF_place_cells_shifted, loc_shifted, BitVector(train_on .& include_mask_shifted), n_pos, true)
            
            
            var, pred, valid = Decoder.direct_basis_decoder(density_map, neuron_map, A_dF_place_cells_binned[i:i+decoding_window, :], n_pos)

            rolling_predicted[i:i+decoding_window, :] .= pred
            rolling_var[i:i+decoding_window] .= var
            rolling_valid[i:i+decoding_window] .= valid
        end
        
        return vcat(rolling_predicted, fill(NaN, shift, 2)), vcat(rolling_var, fill(NaN, shift)), vcat(rolling_valid, fill(NaN, shift))
    end




function quick_decoder(include_mask, A_dF_place_cells, loc_digital, bins::Int, shift::Int, n_pos::Int; decoding_window=120)

        n_sweeps = size(A_dF_place_cells, 1)

        rolling_predicted = fill(NaN, n_sweeps-shift, 2)
        rolling_valid = fill(NaN, n_sweeps-shift)
        rolling_var = fill(NaN, n_sweeps-shift)


        #shift activity against location
        A_dF_place_cells_shifted, loc_shifted, include_mask_shifted = Decoder.shift_time(shift, A_dF_place_cells, [loc_digital, include_mask])

        for_density_maps = copy(Decoder.make_basis_maps(A_dF_place_cells_shifted, loc_shifted, include_mask_shifted, n_pos, true))
        
        density_map = Float32.(sum([for_density_maps[:,:,i] .> nanpctile(for_density_maps[:,:,i], 95) for i in 1:size(for_density_maps, 3)]))
        density_map[density_map .== 0] .= 1

        density_map .= density_map .^(1/4)
        
    
        A_dF_place_cells_binned = Decoder.bin_activity(bins, A_dF_place_cells_shifted)
        
        train_on = trues(n_sweeps-shift)
        index_min = floor(Int32, n_sweeps*3/10)
        index_max = floor(Int32, n_sweeps*6/10)
    
        train_on[index_min:index_max] .= false
        neuron_map = Decoder.make_basis_maps(A_dF_place_cells_shifted, loc_shifted, BitVector(train_on .& include_mask_shifted), n_pos, true)


        var, pred, valid = Decoder.direct_basis_decoder(density_map, neuron_map, A_dF_place_cells_binned[.!train_on, :], n_pos)

        rolling_predicted[index_min:index_max, :] .= pred
        rolling_var[index_min:index_max] .= var
        rolling_valid[index_min:index_max] .= valid
        
        return vcat(rolling_predicted, fill(NaN, shift, 2)), vcat(rolling_var, fill(NaN, shift)), vcat(rolling_valid, fill(NaN, shift))
    end

    function prepare_masks(exclude, n_sweeps, valid_moving_indices, train_on_minutes=30, seed=1)

        # prepare time mask for training
        test_on = trues(n_sweeps)

        minutes_bins = exclude:120:n_sweeps
        minutes = np.digitize(collect(exclude:n_sweeps), minutes_bins)
        Random.seed!(seed)
        train_on = Distributions.sample(1:maximum(minutes)-1, train_on_minutes, replace=false)

        test_on[exclude .+ vcat([findall(minutes .== i) for i in train_on]...)] .= false

        include_mask = trues(n_sweeps)
        include_mask[1:exclude] .= false
        include_mask .= include_mask .& valid_moving_indices
    
        predict_mask = test_on .& include_mask
        train_mask = .!test_on .& include_mask
        
        return predict_mask, train_mask, include_mask

    end


    function shift_time(shift, A_dF, rest)

    rest = [x[1:end-shift] for x in rest]

        return A_dF[shift+1:end, :], rest...
    end


    function bin_activity(bins, A_place_cells)
    
        A_place_cells_binned = Float32.(fill(NaN, size(A_place_cells)))
        for i in 1:size(A_place_cells, 2)
            for j in 1+bins:size(A_place_cells, 1)-bins
                A_place_cells_binned[j, i] = nanmean(A_place_cells[j-bins:j+bins, i])
            end
        end
        return A_place_cells_binned
    end


    function make_basis_maps(A_dF, loc, time_mask, n_pos, filter=false)
    #creates maps for all given activities

        neuron_map_all = fill(NaN32, n_pos, n_pos, size(A_dF,2))

        for neuron_idx in 1:size(A_dF, 2)

            neural_activity, which_loc = MAP.valid_activity_loc(A_dF[:, neuron_idx], time_mask, loc)

            neuron_map_all[:,:,neuron_idx], counted, summed = MAP.calculate_map_direct(neural_activity, which_loc, n_pos; at_least_visit = 1, use_gaussian_filter=filter, sigma=1, filter_mask = nothing)
            neuron_map_all[:,:,neuron_idx] = (neuron_map_all[:,:,neuron_idx] .- nanmean(neuron_map_all[:,:,neuron_idx])) ./ nanstd(neuron_map_all[:,:,neuron_idx])
        end

        return neuron_map_all
    end


    function direct_basis_decoder(density, maps, A_dF_place_cell, n_pos)

        n_sweeps = size(A_dF_place_cell, 1)
        n_neurons = size(A_dF_place_cell, 2)
    
        predicted_pos = fill(NaN, n_sweeps, 2)
        variance = fill(NaN, n_sweeps)
    
        #pre-allocate
        temp = fill(0.0, n_pos, n_pos, 2)
        temp2 = fill(0.0, n_pos, n_pos, 2)
        m = fill(0.0, n_pos, n_pos, 1)
    
        for time in 1:n_sweeps #predict all timepoints - filter out the training and invalid later

            temp .= 0
            temp2 .= 0
        

            threshold = nanpctile(A_dF_place_cell[time, :], 70) #was 70 before -> maybe remove?

            for n in 1:n_neurons
                if A_dF_place_cell[time, n] < threshold
                    continue
                end
                temp[:, :, 2] .= (A_dF_place_cell[time, n] .* (maps[:, :, n]))
                temp[:, :, 1] .= nansum(temp, dim=3)
            end

            m .= copy(temp[:,:,1] ./ density)
            
            m[m .< nanpctile(m, 99)] .= NaN
            variance[time], predicted_pos[time, 1], predicted_pos[time, 2] = spatial_variance(m)
            
        end


        valid = nansum(collect(.!isnan.(A_dF_place_cell[:, :])), dims=2)

        return variance, predicted_pos, valid
    end


    function get_top_neurons(n_place_neurons, population_z, specificity_z)
        # threshold defined by number of neurons
        n = copy(n_place_neurons);
        place_cell_neurons = []
        
    
        mask = .! (isnan.(specificity_z) .| isnan.(population_z))
        
        if length(population_z[mask]) < n_place_neurons
            return findall(mask)
        end

        while length(place_cell_neurons) < n_place_neurons

            top_specificity_z = sortperm(specificity_z[mask],rev=true)[1:n];
            top_population_z = sortperm(population_z[mask],rev=true)[1:n];

            place_cell_neurons = intersect(top_specificity_z,top_population_z)
            n += 1;
        end
        place_cell_neurons = findall(mask)[place_cell_neurons]
    
        if sum(isnan.(population_z[place_cell_neurons])) > 1
            return []
        end
    
        return place_cell_neurons
    end


    function spatial_variance(map) # you might want to filter out low values first
        """
        map: 2d matrix, activity map, mean activity in each bin

        returns: float32, spatial variance
                 float32, x position of centroid
                 float32, y position of centroid

        """
        test_map = copy(map)

        test_map_xy = findall(isfinite.(test_map)); 

        test_map_x = [test_map_xy[i][1] for i = 1:length(test_map_xy)];
        test_map_y = [test_map_xy[i][2] for i = 1:length(test_map_xy)];

        test_map_p = test_map ./ nansum(test_map);
        test_map_p = test_map_p[test_map_xy];


        test_map_mu_x = sum(test_map_x .* test_map_p)
        test_map_mu_y = sum(test_map_y .* test_map_p)

        map_var = sum(((test_map_x .- test_map_mu_x).^2 + (test_map_y .- test_map_mu_y).^2) .* test_map_p)
        return map_var, test_map_mu_x, test_map_mu_y
    end


    function baseline_error(x, y, mask; sample=1000) #x and y already in bin units
        # baseline for random prediction
        trials = fill(NaN, sample)
        
        for i in 1:sample
        
            pred_x = rand(Uniform(minimum(x), maximum(x)), length(x))
            pred_y = rand(Uniform(minimum(y), maximum(y)), length(y))
            
            trials[i] = nanmedian(get_distance(pred_x, x, pred_y, y)[mask])
        end
        return mean(trials), std(trials)
    end


    #returns the x and y of the max value bin
    function get_max_bin(place_map)
        x_max = NaN
        y_max = NaN
        max=-Inf
        for i in 1:size(place_map, 1)
            for j in 1:size(place_map, 2)

                if (place_map[i,j] > max)
                    x_max=i
                    y_max=j
                    max=place_map[i,j]
                end
            end
        end
            return x_max, y_max
    end

    function sort_maps(A_dF_place_cells, loc_digital, n_pos, mask_valid, shift=4)

        A_dF_place_cells_shifted, loc_shifted, mask_shifted = Decoder.shift_time(shift, A_dF_place_cells, [loc_digital, mask_valid])

        place_maps = Decoder.make_basis_maps(A_dF_place_cells_shifted, loc_shifted, mask_shifted, n_pos, true)

        peak_loc_map =  fill(NaN32,2, size(A_dF_place_cells, 2))

        for i in 1:size(A_dF_place_cells, 2)
            peak_loc_map[:,i] .= find_com(place_maps[:,:,i])
        end

        p_x = sortperm(peak_loc_map[1,:])
        p_y = sortperm(peak_loc_map[2,:])

        return p_x, p_y
    end


    function cartesian_product(x,y)
        leny=length(y)
        lenx=length(x)
        m=leny*lenx
        OUT = zeros(Int32, m,2)
        c=1
        for i = 1:lenx
            for j = 1:leny
                OUT[c,1] = x[i]
                OUT[c,2] = y[j]
                c+=1
            end
        end 
        return OUT
    end

    #returns the mean x and y value of multiple top ranked bins
    function get_max_bins(place_map, amount=5)
        x_max = fill(0, amount)
        y_max = fill(0, amount)

        for i in 1:amount
            x_max[i], y_max[i] = Int.(get_max_bin(place_map))
            place_map[x_max[i], y_max[i]] = -Inf
        end

        return mean(x_max), mean(y_max)
    end


    # a way to do useful modulo in Julia, that works with arrays (jumps straight to 1 when too big and to the max when <1)
    function cool_mod(x, m)
        while x < 1
            x = x+m
        end
        while x > m
            x = x-m
        end
        return x
    end

    function bin_to_mm(x)
        return x .* (50/15)
    end


    function get_mean_deviation(x, y)
        return abs.(x .- y)
    end


    function get_distance(x1, y1, x2, y2)
        return sqrt.(abs.(x1 .- y1).^2 + abs.(x2 .- y2).^2)
    end


    function check_cells_component(maps)
        components = fill(NaN, size(maps, 3))
        for i in 1:size(maps, 3)
            components[i] = find_components(maps[:,:,i])
        end
        return components
    end


    function find_components(map)
        components = fill(NaN, 100)

        for p in 1:100
            test_map = copy(map)
            test_map[test_map.<nanpctile(reshape(test_map, length(test_map)), p)].=NaN
            test_map[.!isnan.(test_map)].=1
            test_map[isnan.(test_map)].=0
            img_label = label_components(Int64.(test_map));
            components[p] = nanmaximum(img_label)
        end
        return round(Int32, mean(components[60:90]))
    end
end