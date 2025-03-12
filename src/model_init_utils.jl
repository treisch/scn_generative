module temp_vat_initialize_run

using Distributions, StatsBase, CSV

export load_init_pnw


function load_init_pnw(path_nl, path_el)
    
    # read nodelist:

    nodelist = Dict{Int64, Vector{Int64}}()
    nace_to_consecutive = Dict()
    for row in CSV.File(path_nl)
        push!(nodelist, (row[:firm_id] => [row[:nace_consecutive]+1,0,0,0])) # +1 bc of 0-indexing in pytohn!
        push!(nace_to_consecutive, (row[:nace] => row[:nace_consecutive]+1))
    end

    nace_consecutive_original_counts = Dict((nace_cons=>0) for (nace, nace_cons) in nace_to_consecutive)

    for (k,v) in nodelist
        nace_consecutive_original_counts[v[1]] += 1
    end

    Nnace = maximum(values(nace_to_consecutive))

    # read edgelist
    edgelist = Set{Tuple{Int64, Int64}}()
    for row in CSV.File(path_el)
        push!(edgelist, (row[:supplier], row[:buyer]))
    end

    # initial degree counting
    for k in edgelist
        i,j = k
        nodelist[i][[3,4]] .+= 1  # kout, ktot
        nodelist[j][[2,4]] .+= 1  # kin, ktot
    end

    # remove nodes w/o links
    for k in keys(nodelist)
        if nodelist[k][4]==0
            pop!(nodelist, k)
        end
    end

    Nnodes0 = length(nodelist)
    Llinks0 = length(edgelist)
    
    consecutive_to_nace = Dict((c=>n) for (n,c) in nace_to_consecutive);
    
    return nodelist, edgelist, nace_consecutive_original_counts, nace_to_consecutive, consecutive_to_nace, Nnace, Nnodes0, Llinks0
end



function load_att_mat(path_att_mat, Nnace, nace_to_consecutive)
    am2_ = zeros(Nnace, Nnace);

    for row in CSV.File(path_att_mat)
        i_supplier = nace_to_consecutive[Int(row[:partner_nace2])]
        j_buyer = nace_to_consecutive[Int(row[:subject_nace2])]
        counts = row[:counts]
        am2_[j_buyer, i_supplier] = counts
    end

    am2_ .+= 1e-8;

    # Plots.heatmap(log.(1 .+ am2))

    # const am2 = rand(Nnace, Nnace)
    return am2_ ./ sum.(eachrow(am2_))
end


function load_att_mat_T(path_att_mat, Nnace, nace_to_consecutive)
    am2_ = zeros(Nnace, Nnace);

    for row in CSV.File(path_att_mat)
        i_supplier = nace_to_consecutive[Int(row[:partner_nace2])]
        j_buyer = nace_to_consecutive[Int(row[:subject_nace2])]
        counts = row[:counts]
        am2_[j_buyer, i_supplier] = counts
    end

    am2_ .+= 1e-8;
    
    am2_ = transpose(am2_)

    # Plots.heatmap(log.(1 .+ am2))

    # const am2 = rand(Nnace, Nnace)
    return am2_ ./ sum.(eachrow(am2_))
end

end