
function ax_to_vec(ax)
    (r,c) = size(ax)
    save = []
    for j in 1:r
        for i in 1:c
            append!(save, [ax[j,i]])
        end
    end
    return save
end

function plot_loc(which_neuron; c="r", cmap="viridis", valid_neurons= valid_neurons, create_new=true, X_all = X_all, Y_all = Y_all, Z_all = Z_all, label="", s=1, ax=nothing)
    if create_new
        fig,ax = subplots(2,2, figsize=(12,8))
    end

    ax[1].axis("off")
    
    ax[3].scatter(X_all[valid_neurons], Y_all[valid_neurons], c="gray", s=1)
    ax[3].scatter(X_all[which_neuron], Y_all[which_neuron], c=c, s=s, cmap=cmap)
    ax[3].set_xlabel("X")
    ax[3].set_ylabel("Y")
    ax[3].set_title(label)
    


    ax[2].scatter(X_all[valid_neurons], -Z_all[valid_neurons], c="gray", s=1)
    ax[2].scatter(X_all[which_neuron], -Z_all[which_neuron], c=c, s=s, cmap=cmap)
    ax[2].set_xlabel("X")
    ax[2].set_ylabel("Z")
    ax[2].set_aspect(2)


    ax[4].scatter(Y_all[valid_neurons], -Z_all[valid_neurons], c="gray", s=1)
    ax[4].scatter(Y_all[which_neuron], -Z_all[which_neuron], c=c, s=s, cmap=cmap)
    ax[4].set_xlabel("Y")
    ax[4].set_ylabel("Z")
    ax[4].set_aspect(2)
end

function plot_quantity(quantity; valid_neurons =valid_neurons, c="r", cmap="viridis", vmin = nothing, vmax=nothing, label=nothing, create_new=true, X_all = X_all, Y_all = Y_all, Z_all = Z_all , plot_abs = false)
    p = sortperm(quantity)
    if plot_abs
        p = sortperm(abs.(quantity))
    end
    quantity_r = quantity[p]
    which_neuron_r = valid_neurons[p]
    
    if create_new
        fig,ax = subplots(2,2, figsize=(12,8))
    end

    ax[1].hist(quantity, bins=100)
    ax[1].set_title(label)
    
    ax[3].scatter(X_all[which_neuron_r], Y_all[which_neuron_r], c=quantity_r, s=2, cmap=cmap, vmin=vmin, vmax= vmax)
    ax[3].set_xlabel("X")
    ax[3].set_ylabel("Y")    


    ax[2].scatter(X_all[which_neuron_r], -Z_all[which_neuron_r],c=quantity_r, s=2, cmap=cmap, vmin=vmin, vmax= vmax)
    ax[2].set_xlabel("X")
    ax[2].set_ylabel("Z")
    ax[2].set_aspect(2)


    ax[4].scatter(Y_all[which_neuron_r], -Z_all[which_neuron_r], c=quantity_r, s=2, cmap=cmap, vmin=vmin, vmax= vmax)
    ax[4].set_xlabel("Y")
    ax[4].set_ylabel("Z")
    ax[4].set_aspect(2)
end

py"""
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
"""


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


function activity_trajectory(neural_activity, x_fish_valid, y_fish_valid,alpha)
    fig = figure(dpi=250)
    plot(x_fish_valid, y_fish_valid, "k.", markersize = 0.1)
    
    layer_rank = sortperm(neural_activity[1:end-1])
    for i = layer_rank
        cur_value = neural_activity[i]; 
        cur_value < 0 ? cur_value = 0 : (cur_value > 1 ? cur_value = 1 : 1+1)
        if isfinite(cur_value)
            plot(x_fish_valid[i:i+1], y_fish_valid[i:i+1], color = [cur_value, 0, 0], alpha = 0.5, linewidth = 2)
        else
            plot(x_fish_valid[i:i+1], y_fish_valid[i:i+1], color = [0, 0, 0], alpha = 0.1)
        end

    end
    imshow(img_bg_end', cmap="binary", origin="lower",alpha=alpha)
    axis("off")
    return fig
end


function plot_mean_density(bins, y_tog; plot_color="k", shade_color="gray", line_color="k", plot_median = true, label= "", zorder = 0)
    
    n_bins = length(bins)
    x = (bins[1:end-1] + bins[2:end])/2
    y = mean(y_tog)
    ymax = y + std(y_tog)
    ymin = y - std(y_tog);
    
    cumsum_y_mean = cumsum(mean(y_tog)*(maximum(bins) - minimum(bins))/(n_bins-1))

    x_median = numpy.round(numpy.interp(0.5, cumsum_y_mean, x), 2)
    x_14 = numpy.round(numpy.interp(0.25, cumsum_y_mean, x), 2)
    x_34 = numpy.round(numpy.interp(0.75, cumsum_y_mean, x), 2)
    
    
    plot(x, y, color= plot_color, label= label, zorder = zorder)
    fill_between(x, ymax, ymin, color=shade_color, alpha=0.3, zorder = zorder)
    if plot_median
        axvline(x_median, color=line_color, label= "$x_median", linestyle="dashed")
        # axvline(x_median, color=line_color, label= "$x_median ($x_14 - $x_34)")
    end
    
    

end