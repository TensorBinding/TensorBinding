# -*- coding: utf-8 -*-
using ITensors: MPO, MPS, OpSum, expect, inner, siteinds
using ITensors
using ITensorMPS
using LinearAlgebra
using Plots
using QuanticsTCI
import TensorCrossInterpolation as TCI
using Quantics
using StatsBase
using Base.Threads

function get_error_plot(conv_error)
 
    # Prepare data
    x = 1:length(conv_error)  # Indices for the x-axis
    y = conv_error            # Error values for the y-axis

    # Plot the data with exponential scaling on the y-axis
    plot(
        x, y,
        yscale = :log10,            # Set y-axis to logarithmic scale (base 10 exponential)
        xlabel = "Iteration number",          # Label for x-axis
        ylabel = "Convergence Error", # Label for y-axis
        title = "Convergence Error Plot", # Title for the plot
        legend = false,            # Disable legend
        lw = 2                     # Line width
    )
end

"""
For random vectors used for stochastic tracing
"""
function slice_edges(L, n)
    step = round(Int, L / n)  # Round L/n to the nearest integer
    edges = [(i * step) for i in 0:(L - 1) ÷ step]  # Calculate edges
    if edges[end] < L - 1
        push!(edges, L - 1)  # Ensure the final edge is L - 1
    end
    return edges
end

function generate_complex_random(num,L,sites)
     
    the_mps = (rand() - 0.5) * exp.(2im * π * rand()) * random_mps(sites, to_binary_vector(Int(num), L))
  
    return  the_mps 
end
# -

function generate_mps_list(L_ture, L, num_frag, sites,r_num,rr)
    # Get the edges of the main interval
    edges = slice_edges(L_ture, num_frag)
    mps_list = []
    
    # Iterate through each interval
    for i in 1:length(edges) - 1
        rand_list = []
        start = edges[i]
        stop = edges[i + 1]
 
        interval_size = stop - start
        num_to_sample = min(rr, interval_size)  # Ensure we don't sample more than the interval size
        num_random_values = sort(StatsBase.sample(start:stop-1, num_to_sample; replace = false))
        for j in 1:r_num
            # Array to store results from each thread
            r_vecs = [0 * random_mps(sites) for _ in 1:nthreads()]

            @threads for num in num_random_values#start:stop-1
                
                r_vecs[threadid()] = +( r_vecs[threadid()], generate_complex_random(Int(num), L,sites);maxdim = 90)
            end

            # Combine all r_vecs from each thread
            r_vec = sum(r_vecs)
            r_vec = r_vec/norm(r_vec)
            r_vec= ITensorMPS.truncate!(r_vec;maxdim = 90)
            push!(rand_list, r_vec)
        end
   
        # Add the inner list to the main list
        push!(mps_list, rand_list)
    end
    GC.gc()
    return mps_list
end

"""
for calculation of Chebyshev moments
"""
function get_mus_combine(tn_folder, mps_folder, sites, L, num_frag, N)
   
    # Step 1: Get sorted file lists
    tn_files = sort(readdir(tn_folder, join=true), by=file -> parse(Int, match(r"(\d+)", basename(file)).match))
    mps_files = sort(readdir(mps_folder, join=true), by=file -> parse(Int, match(r"(\d+)", basename(file)).match))

    # Step 2: Initialize mus list
    mus = []

    # Step 3: Compute mus for each Tn
    for (i, tn_file) in enumerate(tn_files)
        # Load Tn
        Tn = h5open(tn_file, "r") do f
            read(f, "hampo",MPO)  # Assuming Tn is stored as "Tn" in the file
        end
         
        # Initialize mu_list for current Tn
        mu_list = []

        # Compute inner products with all MPS
        for mps_file in mps_files
            # Load MPS
            vec = h5open(mps_file, "r") do f
                read(f, "random_vec", MPS)  # Assuming MPS is stored as "mps" in the file
            end
           # vec= ITensorMPS.truncate!(vec;maxdim = 90)
            # Compute inner product <r|Tn|r> for each vec in vec_list
            mu_r =  inner(vec', Tn, vec) 
            push!(mu_list, mu_r)

            # Delete MPS file after use
            vec = nothing
            GC.gc()  # Force garbage collection
        end

        # Store mu_list for current Tn
        push!(mus, mu_list)

        # Delete Tn file after use
        Tn = nothing
 
        GC.gc()
    
    end

    # Step 4: Apply Jackson kernel
    jackson_kernel = [(N - n) * cos(π * n / N) + sin(π * n / N) / tan(π / N) for n in 0:N-1]
    jackson_kernel /= N

    for i in 1:N
        mus[i] = jackson_kernel[i] .* mus[i]
    end

    return mus
end

#for smaller moire
function get_mus_small(tn_folder,  sites, L, N)

    # Step 1: Get sorted file lists
    tn_files = sort(readdir(tn_folder, join=true), by=file -> parse(Int, match(r"(\d+)", basename(file)).match))
   

    # Step 2: Initialize mus list
    mus = []

    # Step 3: Compute mus for each Tn
    for (i, tn_file) in enumerate(tn_files)
        # Load Tn
        Tn = h5open(tn_file, "r") do f
            read(f, "hampo",MPO)  # Assuming Tn is stored as "Tn" in the file
        end
         
        # Initialize mu_list for current Tn
        mu_list = []

        # Compute inner products with all MPS
        #we pick the middle 100 sites
        for sit_num in 2^(L - 1)-50:2^(L - 1) + 50
            # Load MPS
            vec =  random_mps(sites, to_binary_vector(Int(sit_num), L))
            mu_r =  inner(vec', Tn, vec) 
            push!(mu_list, mu_r)

            # Delete MPS file after use
            vec = nothing
            GC.gc()  # Force garbage collection
        end

        # Store mu_list for current Tn
        push!(mus, mu_list)

        # Delete Tn file after use
        Tn = nothing
        
        println("d")
        GC.gc()
    end

    # Step 4: Apply Jackson kernel
    jackson_kernel = [(N - n) * cos(π * n / N) + sin(π * n / N) / tan(π / N) for n in 0:N-1]
    jackson_kernel /= N

    for i in 1:N
        mus[i] = jackson_kernel[i] .* mus[i]
    end

    return mus
end

function get_ldos(mus_up, mus_down, num_frag,en_num)

    xvals = -num_frag/2:num_frag/2 - 1
     
    yvals = range(-0.99, 0.99; length=en_num)
    ldos_matrix_up = zeros(en_num,num_frag)  
    ldos_matrix_down = zeros(en_num,num_frag)

    # Compute ldos for each y and k, summing over all mus elements
    for k in 1:N
        
        for (i, y) in enumerate(yvals)

            if k == 1
                for m in 1:num_frag
                    ldos_matrix_up[i, m] +=  real.(mus_up[k][m]) /(π* sqrt(1 - y^2 ))
                    ldos_matrix_down[i, m] += real.(mus_down[k][m]) /(π* sqrt(1 - y^2 ))
                end
                
            else  
            # Calculate the scalar coefficient
                coefficient = cos((k-1) * acos(y))  /(π* sqrt(1 - y^2 ))    
                # Sum contributions from all `mus` elements
                for m in 1:num_frag
                    ldos_matrix_up[i, m] += 2 * real.(mus_up[k][m]) * coefficient
                    ldos_matrix_down[i, m] += 2 * real.(mus_down[k][m]) * coefficient
                end
            end
        end
    end
    
    ldos_matrix = ldos_matrix_up + ldos_matrix_down 

    nor = maximum(ldos_matrix) 
    ldos_matrix= ldos_matrix/nor#normalize
    #plot ldos
    
    quantity = floor(log10(2^(L - 1)))

    # Identify indices for xticks near the edges and 0
    edge_offset = Int(round(length(xvals) * 0.1) ) # Choose 10% from the edges as offset
    index_left = Int(round(edge_offset))
    index_right = Int(round(length(xvals) - edge_offset))
    index_center = Int(round(length(xvals) / 2))

    # Map indices to values
    selected_xvals = [xvals[index_left], 0, xvals[index_right]]

    # Generate labels for the xticks
    formatted_labels = [
        x == 0 ? "0" : " $(x < 0 ? "-" : "")5 × 10^{$(Int(quantity))}" for x in selected_xvals
    ]
 

    #A quick showing
    plot(size=(800, 400))
    p = contourf(
        xvals,
        yvals .* 10,
        real.(ldos_matrix),
        xlabel="site",
        ylabel="Energy",
        title="LDOS Contour Plot from MPO",
        levels=100,
        c=:inferno,
        linecolor=:transparent,
        linewidth=0.0,
        xticks=(selected_xvals, formatted_labels),
        tick_direction=:out,
        ylims=(-4, 4),  
        colorbar_title="DOS (normalized)"
        
    )
     
    
    display(p)
    return ldos_matrix
end

function get_dos(ldos_matrix)
    
    dos =   sum(ldos_matrix, dims=2)  
    nor = maximum(dos)
    dos = dos / nor
    xvals = 10* range(-0.99, 0.99; length=length(dos))

    p = plot(
        xvals, dos,
        xlabel = "Energy",          # Label for x-axis
        ylabel = "DOS (normalized)", # Label for y-axis
        title = "DOS Plot from MPO", # Title for the plot
        legend = false,            # Disable legend
        lw = 2                     # Line width
    )
    display(p);
    return dos
end

#not applicable for large system
function get_magnetism(qtt_up, qtt_down, L)
    mz = []
    xvals = 1:(2^L)  
    xv  = -(2^L)/2:(2^L)/2 - 1
    
    for x in xvals 
        z = abs.(qtt_den_new_up(Int(x)) - qtt_den_new_down(Int(x)))
        push!(mz,z)
    end

    up_bound = maximum(mz) + 0.05
    low_bound = minimum(mz) - 0.05
    p = plot(
        xvals , mz,
        xlabel = "Site numebr",          # Label for x-axis
        ylabel = "|Magnetism|", # Label for y-axis
        title = "Magnetism vs site number", # Title for the plot
        legend = false,            # Disable legend
        lw = 2,                     # Line width
        ylim = (low_bound,up_bound )
    )
    display(p)
    return mz
end
