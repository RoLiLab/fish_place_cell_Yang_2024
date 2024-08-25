

function example_stuff()
    A_dFF_place_cells = A_dFF[:, place_candidates_unique]
    

    rolling_predicted, rolling_var, rolling_valid = Decoder.rolling_decoder(bool_index, A_dFF_place_cells, loc_digital, activity_bins, activity_shift, n_pos)
    
    errors = Decoder.get_distance(rolling_predicted[:, 1], x_in_bins, rolling_predicted[:, 2], y_in_bins)
    print(nanmedian(errors[bool_index]).*long_axis_in_mm/long_axis_in_bins, " mm")
    
    
    sorting_x, sorting_y = Decoder.sort_maps(A_dFF_place_cells, loc_digital, n_pos, bool_index)
    
    h5open(joinpath(data_path(ds_Lorenz), "$(experiment_filename)$(file_name)"), "w") do file
        file["place_cell_included"] = place_candidates_unique
        file["x_in_bins"] = x_in_bins
        file["y_in_bins"] = y_in_bins
        file["x_predicted_in_bins"] = rolling_predicted[:, 1]
        file["y_predicted_in_bins"] = rolling_predicted[:, 2]
        file["bool_index"] = Int8.(bool_index)
        file["moving_valid"] = Int8.(moving_valid)
        file["activity_sorting_x"] = sorting_x
        file["activity_sorting_y"] = sorting_y
        file["img_bg_end"] = img_bg_end
    end
end



function error_by_populations(use_neurons, populations)

    errors_mm = Dict()
    for name in keys(populations)

        A_dFF_place_cells = A_dFF[:, populations[name]]
        rolling_predicted, rolling_var, rolling_valid = Decoder.rolling_decoder(bool_index, A_dFF_place_cells, loc_digital, activity_bins, activity_shift, n_pos)

        errors = Decoder.get_distance(rolling_predicted[:, 1], x_in_bins, rolling_predicted[:, 2], y_in_bins)

        errors_mm[name] = nanmedian(errors[bool_index]) * long_axis_in_mm/long_axis_in_bins
    end
    
    
    
    fig, ax = plt.subplots()
    ax.barh(y=1:length(errors_mm), width=values(errors_mm))
    ax.set_yticks(1:length(errors_mm), labels=keys(errors_mm))
    ax.set_xlabel("error [mm]")
    
    
    baseline_mean, baseline_std = Decoder.baseline_error(x_digital, y_digital, bool_index; sample=use_neurons)
    errors_mm["baseline_uniformly_random"] = baseline_mean * long_axis_in_mm/long_axis_in_bins
    
    
    neural_activity, which_loc = MAP.valid_activity_loc(A_dFF[:, 1], bool_index, loc_digital)
    neuron_map, counted, summed = MAP.calculate_map_direct(neural_activity, which_loc, n_pos; at_least_visit = 1, use_gaussian_filter=false, sigma=1, filter_mask = nothing)
    var, com_x, com_y = Decoder.spatial_variance(counted)
    errors_mm["baseline_behavior_informed"] = nanmedian(Decoder.get_distance(x_in_bins, fill(com_x, length(x_in_bins)), y_in_bins, fill(com_y, length(y_in_bins)))[bool_index]) * long_axis_in_mm/long_axis_in_bins
    
    
    h5open(joinpath(data_path(ds_Lorenz), "$(experiment_filename)$(file_name)"), "cw") do file
        try
            file["errors_by_population_in_mm_keys"] = Vector{String}(collect(keys(errors_mm)))
            file["errors_by_population_in_mm_values"] = Vector{Float32}(collect(values(errors_mm)))
        catch
            delete_object(file, "errors_by_population_in_mm_keys")
            delete_object(file, "errors_by_population_in_mm_values")
            file["errors_by_population_in_mm_keys"] = Vector{String}(collect(keys(errors_mm)))
            file["errors_by_population_in_mm_values"] = Vector{Float32}(collect(values(errors_mm)))
        end
    end
end

    
    
    
    
function error_by_amount(;space_n=120) #based on specificity

    amount = [10000, 5000, 1000, 500, 250, 100, 50, 25, 10, 5, 4, 3, 2, 1] 

    #forebrain
    sample_cells = [mask_tel[Decoder.get_top_neurons(x, specificity_population_z[mask_tel], specificity_shuffle_z[mask_tel])] for x in amount]


    error_by_amount_place_cells = fill(NaN, length(sample_cells))

    @showprogress for i in 1:length(sample_cells)

        if sample_cells[i] == []
            error_by_amount_place_cells[i] = NaN
            continue
        end

        A_dFF_place_cells = copy(A_dFF[:, sample_cells[i]])

        rolling_predicted, rolling_var, rolling_valid = Decoder.rolling_decoder(bool_index, A_dFF_place_cells, loc_digital, activity_bins, activity_shift, n_pos; space=space_n)

        errors = Decoder.get_distance(rolling_predicted[:, 1], x_in_bins, rolling_predicted[:, 2], y_in_bins)
        error_by_amount_place_cells[i] = nanmedian(errors[bool_index])

    end
    
    plt.plot(amount, error_by_amount_place_cells.*(long_axis_in_mm/long_axis_in_bins))
    plt.xlabel("cells included")
    plt.ylabel("error [mm]")
    plt.xscale("log")


    h5open(joinpath(data_path(ds_Lorenz), "$(experiment_filename)$(file_name)"), "cw") do file
        try
            file["errors_by_amount_$(space_n)_in_mm_keys"] = amount
            file["errors_by_amount_$(space_n)_in_mm_values"] = error_by_amount_place_cells.*(long_axis_in_mm/long_axis_in_bins)
        catch
            delete_object(file, "errors_by_amount_$(space_n)_in_mm_keys")
            delete_object(file, "errors_by_amount_$(space_n)_in_mm_values")
            file["errors_by_amount_$(space_n)_in_mm_keys"] = amount
            file["errors_by_amount_$(space_n)_in_mm_values"] = error_by_amount_place_cells.*(long_axis_in_mm/long_axis_in_bins)
        end
    end
end


function error_by_amount_greedy(amount=25) # based on finding the best next cell to include in decoding

    # get the neurons
    neurons_chosen = Int[]

    left_neurons = copy(place_candidates_unique)

    errors = Float32[]

    @showprogress for i in 1:amount # include more neurons

        current_min = Inf
        current_argmin = NaN

        @showprogress for n in left_neurons
            A_dFF_candidate_cells = @view A_dFF[:, cat(neurons_chosen, [n], dims=1)]

            # exact version to get the exact result
            #pred2, var2, valid2 = Decoder.rolling_decoder(bool_index, A_dFF_candidate_cells, loc_digital, activity_bins, activity_shift, n_pos; decoding_window=120);

            # fast version to get a quick result
            pred2, var2, valid2 = Decoder.quick_decoder(bool_index, A_dFF_candidate_cells, loc_digital, activity_bins, activity_shift, n_pos; decoding_window=120)

            errors2 = Decoder.get_distance(pred2[:, 1], x_in_bins, pred2[:, 2], y_in_bins);

            if nanmedian(errors2[bool_index]) < current_min
                current_argmin = copy(n)
                current_min = copy(nanmedian(errors2[bool_index]))
            end
        end

        append!(errors, [current_min])
        append!(neurons_chosen, [current_argmin])
    end
    
    
    # do the proper decoding based on the cells found
    errors = []
    @showprogress for i in 1:length(neurons_chosen) # include more neurons

        current_min = Inf
        current_argmin = NaN

        A_dFF_candidate_cells = @view A_dFF[:, neurons_chosen[1:i]]

        pred2, var2, valid2 = Decoder.rolling_decoder(bool_index, A_dFF_candidate_cells, loc_digital, activity_bins, activity_shift, n_pos; decoding_window=120);

        errors2 = Decoder.get_distance(pred2[:, 1], x_in_bins, pred2[:, 2], y_in_bins);


        append!(errors, [nanmedian(errors2[bool_index])])
    end
    
    
    plt.plot(1:length(neurons_chosen), errors.*(long_axis_in_mm/long_axis_in_bins))
    plt.xlabel("cells included")
    plt.ylabel("error [mm]")
    title("greedy chosing")
    
    
    
    h5open(joinpath(data_path(ds_Lorenz), "$(experiment_filename)$(file_name)"), "cw") do file
        try
            file["errors_by_amount_greedy_in_mm_values"] = errors.*(long_axis_in_mm/long_axis_in_bins)
            file["errors_by_amount_greedy_in_mm_cells"] =  Int.(neurons_chosen)
        catch
            delete_object(file, "errors_by_amount_greedy_in_mm_values")
            delete_object(file, "errors_by_amount_greedy_in_mm_cells")
            
            file["errors_by_amount_greedy_in_mm_values"] = errors.*(long_axis_in_mm/long_axis_in_bins)
            file["errors_by_amount_greedy_in_mm_cells"] =  Int.(neurons_chosen)
        end
    end
end


function error_bin_shift()

    A_dF_place_cells = A_dFF[:, place_candidates_unique]
    

    bins_shift_mm = fill(NaN, 10,15)

    prod = Decoder.cartesian_product(collect(bins[1]:bins[2]), collect(shifts[1]:shifts[2]))

    @showprogress for p in 1:size(prod, 1)
        bin = prod[p, 1]
        shift = prod[p, 2]
            #shift activity against location

            predicted, rolling_var, rolling_valid = Decoder.rolling_decoder(bool_index, A_dF_place_cells, loc_digital, bin, shift, n_pos)


            errors = Decoder.get_distance(predicted[:, 1] , x_in_bins, predicted[:, 2], y_in_bins)
            bins_shift_mm[bin+1, shift+1] = nanmedian(errors[bool_index])*long_axis_in_mm/n_pos
    end
    
    imshow(bins_shift_mm', origin="lower")

    h5open(joinpath(data_path(ds_Lorenz), "$(experiment_filename)$(file_name)"), "cw") do file
        try
            file["error_by_bins_shift_mm"] = bins_shift_mm
        catch
            delete_object(file, "error_by_bins_shift_mm")
            file["error_by_bins_shift_mm"] = bins_shift_mm
        end
    end
end
    
    
function error_long_shift()
    amount = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 30, 50, 60, 80, 100, 300, 500, 600, 700, 800, 900, 1000]
    long_shift_errors_mm = fill(NaN, length(amount))

    @showprogress for i in 1:length(amount)
        shift = amount[i]
        
        predicted, rolling_var, rolling_valid = Decoder.rolling_decoder(bool_index, A_dFF[:, place_candidates_unique], loc_digital, activity_bins, shift, n_pos)
        
        errors = Decoder.get_distance(predicted[:, 1] , x_in_bins, predicted[:, 2], y_in_bins)

        long_shift_errors_mm[i] = nanmedian(errors[bool_index])*long_axis_in_mm/long_axis_in_bins
    end

    plt.plot(amount, long_shift_errors_mm)
    plt.axhline(20, alpha=0.5)
    plt.xlabel("shift")
    plt.ylabel("error [mm]")
    plt.xscale("log")


    h5open(joinpath(data_path(ds_Lorenz), "$(experiment_filename)$(file_name)"), "cw") do file
        try
            file["long_shift_errors_mm_keys"] = amount
            file["long_shift_errors_mm_values"] = long_shift_errors_mm
        catch
            delete_object(file, "long_shift_errors_mm_keys")
            delete_object(file, "long_shift_errors_mm_values")
            file["long_shift_errors_mm_keys"] = amount
            file["long_shift_errors_mm_values"] = long_shift_errors_mm
        end
    end
end



function decoding_animation(bins, shift)
    #animate decoder performance

    A_dFF_place_cells = A_dFF[:, place_candidates_unique]

    n_sweeps = size(A_dFF_place_cells, 1)


    #shift activity against location
    A_dFF_place_cells_shifted, loc_shifted, include_mask_shifted = Decoder.shift_time(shift, A_dFF_place_cells, [loc_digital, bool_index])

    A_dFF_place_cells_binned = Decoder.bin_activity(bins, A_dFF_place_cells_shifted)

    decoder_maps = fill(NaN, n_pos, n_pos, n_sweeps)

    @showprogress for i in 120:120:n_sweeps-300
        train_on = trues(n_sweeps-shift)
        train_on[i-120+1:i+240] .= false

        neuron_map = Decoder.make_basis_maps(A_dFF_place_cells_shifted, loc_shifted, train_on, n_pos, true)

        for time in i:i+120

            temp=fill(0.0, n_pos, n_pos, 2)


            for n in 1:size(A_dFF_place_cells_binned, 2)
                temp[:, :, 2] = (A_dFF_place_cells_binned[time, n] .* (neuron_map[:, :, n]))
                temp[:, :, 1] = nansum(temp, dim=3)
            end
            decoder_maps[:, :, time] = copy(temp[:, :, 1])# ./ nanstd(temp[:, :, 1])
        end

    end


    h5open(joinpath(data_path(ds_Lorenz), "$(experiment_filename)$(file_name)"), "cw") do file
        try
            file["decoder_maps"] = decoder_maps;
        catch
            delete_object(file, "decoder_maps")
            file["decoder_maps"] = decoder_maps;
        end
    end
end
