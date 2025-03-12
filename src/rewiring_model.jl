module temp_vat
using Random, Distributions, LinearAlgebra, StatsBase, CSV, DataFrames

export run_sim, calc_degree_hist, save_run

function dropLink!(supplier_id, customer_id, nl, el)
    if (supplier_id, customer_id) in el
        pop!(el, (supplier_id, customer_id))
        nl[customer_id][[2,4]] .-= 1 # reduce kin and ktot
        nl[supplier_id][[3,4]] .-= 1 # reduce kout and ktot
    else
        throw(ArgumentError("Link ($supplier_id, $customer_id) doesn't exist."))
    end
    nothing
end

function addLink!(supplier_id, customer_id, nl, el)
    if (supplier_id, customer_id) in el
        throw(ArgumentError("Link ($supplier_id, $customer_id) exists already."))
    else
        push!(el, (supplier_id, customer_id))
        nl[customer_id][[2,4]] .+= 1 # increase kin and ktot
        nl[supplier_id][[3,4]] .+= 1 # increase kout and ktot
    end
    nothing
end

# generate attachment potential vec ktilde
function attachment_potential(dangling_il, nl, fda_)
    return Dict(v[4]!=0 ? k=>fda_(v[4]) :  k=>fda_(dangling_il[k]) for (k,v) in nl) # ktilde
end

function increase_dictcnt!(dct, key, inc)
    if key in keys(dct)
        dct[key] += inc
    else
        push!(dct, (key => inc))
    end
end

# remove links (remember which)
function remove_links_leave_some_dangling!(pterm, prepl, nodelist, edgelist, dangling_inlinks_)
    links = collect(edgelist)
    nremove = Int(floor(length(links)*pterm))
    remove_links = shuffle(links)[1:nremove]

    # remove nodes (remember which)
    for l in remove_links
        dropLink!(l[1], l[2], nodelist, edgelist)
    end

    # replace "dangling" links" (with prep)
    nreplace = Int(floor(nremove*prepl))
    # add to dangling links thing, fix all of them together later.
    if nreplace>0
        for l in shuffle(remove_links)[1:nreplace]
            increase_dictcnt!(dangling_inlinks_, l[2], 1)
        end
    end
end

function remove_nodes_wo_links!(nl, dangling_il, inlink_targets_nn)
    cnt = 0
    for (k,v) in nl
        if (v[4]==0) & !(k in keys(dangling_il)) & !(k in keys(inlink_targets_nn))
            pop!(nl,k)
            cnt+=1
        end
    end
    println("Removed ", cnt, " nodes.")
end

# spontaneously add links
## generate number of new links for every node and add to dangling links
function spontaneously_add_dangling_links!(nl, dangling_il, nsp_, ktm1)
    N_ = length(nl)
    for i in keys(nl)
        ktot_i = ktm1[i]
        nace_i = nl[i][1]
        nsp_i = nsp_(ktot_i, nace_i, N_)[1]
        if nsp_i>0
            increase_dictcnt!(dangling_il, i, nsp_i)
        end
    end
end

function to_supp_view(nodelist, edgelist)
    # creates dict <customer_i> => Set(<all suppliers_j of customer_i>)
    supp_view = Dict((k=>Set{Int64}()) for k in keys(nodelist))
    for l in edgelist
        push!(supp_view[l[2]], l[1])
    end
    return supp_view
end

function to_cust_view(nodelist, edgelist)
    # creates dict <customer_i> => Set(<all suppliers_j of customer_i>)
    cust_view = Dict((k=>Set{Int64}()) for k in keys(nodelist))
    for l in edgelist
        push!(cust_view[l[1]], l[2])
    end
    return cust_view
end

function to_sec_view(nl,Nnace)
    sec_view = Dict((k=>Set{Int64}()) for k in 1:Nnace)
    for (k,v) in nl
        push!(sec_view[v[1]], k)
    end
    sec_view = Dict((k=>collect(v)) for (k,v) in sec_view)
    return sec_view
end

function subset_attachment_potential(ap_, sec_lst)
    Dict((k=>v) for (k,v) in ap_ if k in sec_lst)
end

function exclude_self_and_current_supps!(ap_temp, i_self, supp_view)
    for j in supp_view[i_self]
        if j in keys(ap_temp)
            ap_temp[j]=0
        end
    end
    if i_self in keys(ap_temp)
        ap_temp[i_self]=0
    end
    nothing
end


function sample_new_target(candidate_nodes, ap_temp, ntargets)
    prob_weights = getindex.(Ref(ap_temp), candidate_nodes)
    return sample(candidate_nodes, Weights(prob_weights), ntargets, replace=false)
end

function get_new_target(i_customer, sec_view_, ap_sec_view_, supp_view_, ntargets)
    # candidate nodes must be a list!
#         candidate_nodes = collect(candidate_nodes) # collect(sec_view[2])
    
    ap_temp = deepcopy(ap_sec_view_)
    #
    exclude_self_and_current_supps!(ap_temp, i_customer, supp_view_)
    
    n_candidates = sum([1 for v in values(ap_temp) if v>0])
    if (n_candidates==0)
        println("No egligble targets.")
        return Set{Int64}()
    elseif n_candidates<ntargets
        println("Fewer available targets than needed, reduce ntargets. i=$i_customer, ntargets=$ntargets")
        ntargets = n_candidates
    end

    new_targets = sample_new_target(sec_view_, ap_temp, ntargets)
    return new_targets
end


function connect_dangling_bonds(i_cust, nace_i, n_dangling, Nnace, att_mat_2, sec_view, ap_sec_view, supp_view)
    # pick sector
#     println(Nnace," ", n_dangling," ", nace_i," ", length(att_mat_2[nace_i,:]))
    secvec = countmap(sample(1:Nnace, Weights(att_mat_2[nace_i,:]), n_dangling))
    
    # pick nodes
    targets = Set{Int64}()
    for (sec, n_targets) in secvec
        union!(targets,get_new_target(i_cust, sec_view[sec], ap_sec_view[sec], supp_view, n_targets))
    end
    return targets
end

function remove_nodes_spontaneously!(nl, el, dangling_il, pexit, prepl)
    nremove = Int(floor(length(nl)*pexit))
    remove_nodes = shuffle(collect(keys(nl)))[1:nremove]
    
    for l in el
        if l[1] in remove_nodes
            dropLink!(l[1], l[2], nl, el)
        elseif l[2] in remove_nodes
            if rand()<prepl
                increase_dictcnt!(dangling_il, l[2], 1)
            end
            dropLink!(l[1], l[2], nl, el)
        end
    end
end


function calc_degree_hist(nl)
    kin = Dict()
    kout = Dict()
    ktot = Dict()
    for (k,v) in nl
        increase_dictcnt!(kin, v[2], 1)
        increase_dictcnt!(kout, v[3], 1)
        increase_dictcnt!(ktot, v[4], 1)
    end
    return kin, kout, ktot
end

function Lt(nodelist)
    sum([v[2] for (k,v) in nodelist])
end

function dangling_inlinks_to_sec_view(dangling_inlinks_, nl_, nsec_)
    # how does dangling_il look?
    # (idx => n_dangling)   
    dangling_il_sec_view = Dict{Int64, Dict{Int64,Int64}}((i=>Dict{Int64, Int64}()) for i in 1:nsec_)
    
    for (i, n_dangling) in dangling_inlinks_
        nace_i = nl_[i][1]
        push!(dangling_il_sec_view[nace_i], (i=>n_dangling))
    end
    
    return dangling_il_sec_view
end

function add_outlinks_to_new_nodes!(nl, el, new_nodes_, dangling_inlinks_, am2_T, nsec_, cust_view)
    # prepare: dangling_inlinks sector view!
    dangling_il_sec_view = dangling_inlinks_to_sec_view(dangling_inlinks_, nl, nsec_)
    
    sec_view_ids = Dict((sec=>collect(keys(v))) for (sec,v) in dangling_il_sec_view)
    
    # for every new node
    for (node_id, node_sector_kin_kout)  in new_nodes_
        node_sec = node_sector_kin_kout[1]
        node_kout = node_sector_kin_kout[3]
        
        # sample sector using am2_T
        secvec = countmap(sample(1:nsec_, Weights(am2_T[node_sector_kin_kout[1],:]), node_kout))
    
        # pick nodes
        targets = Set{Int64}()
        for (sec, n_targets) in secvec
            union!(targets, get_new_target(node_id, sec_view_ids[sec], dangling_il_sec_view[sec], cust_view, n_targets))
        end
        
        for j in targets
            if !((node_id, j) in el)
                addLink!(node_id, j,nl,el)
                dangling_inlinks_[j]-=1
#             else
#                 println("Link ($j, $i) already in edgelist_.")
            end
        end
    end
    nothing
end
    
        # sample target from dangling inlinks in that sector
    
        # add link to edgelst, update nodelist, decrease dangling_inlinks_ and d_il_sec_view_
        
        




# timestep function:
function timestep!(nodelist_, edgelist_, pterm, prepl, p_node_exit, nsp, Nnace, am2, am2_T, fda_, new_nodes_fct_)
    idx_it = maximum(collect(keys(nodelist_)))+1
    
    # add nodes from distribution
    new_nodes = new_nodes_fct_(nodelist_)
    inlink_target_new_nodes = Dict()
    new_nodes_dict = Dict{Int64, Vector{Int64}}()
    
    for nn in new_nodes
        push!(nodelist_, (idx_it => [nn[1],0,0,0]))
        push!(inlink_target_new_nodes,(idx_it => nn[2]))
        push!(new_nodes_dict, (idx_it => nn))
        idx_it = idx_it + 1
    end
       
    
    println("Added $(length(new_nodes)) new nodes.")
    
    ktot_tm1 = Dict(v[4]!=0 ? k=>v[4] :  k=>new_nodes_dict[k][2]+new_nodes_dict[k][3] for (k,v) in nodelist_) # ktilde

    # make attachment potential vector
    ap = attachment_potential(inlink_target_new_nodes, nodelist_, fda_)
    
    dangling_inlinks = inlink_target_new_nodes
    
    # remove links
    remove_links_leave_some_dangling!(pterm, prepl, nodelist_, edgelist_, dangling_inlinks)
    
    # remove nodes without links, later also spontaneously remove nodes
    remove_nodes_wo_links!(nodelist_, dangling_inlinks, inlink_target_new_nodes)
    
    # remove node spont.
    remove_nodes_spontaneously!(nodelist_, edgelist_, dangling_inlinks, p_node_exit, prepl)
    remove_nodes_wo_links!(nodelist_, dangling_inlinks, inlink_target_new_nodes)
    
    # spontaneously add dangling links!
    spontaneously_add_dangling_links!(nodelist_, dangling_inlinks, nsp, ktot_tm1)

    supp_view = to_supp_view(nodelist_, edgelist_);
    cust_view = to_cust_view(nodelist_, edgelist_);
    
    sec_view = to_sec_view(nodelist_, Nnace)
    
    ntot_dangling = sum(collect(values(dangling_inlinks)))
    
    # do outlinking of new nodes first:
    add_outlinks_to_new_nodes!(nodelist_, edgelist_, new_nodes_dict, dangling_inlinks, am2_T, Nnace, cust_view)
    
    println("ntot_dangling before: $ntot_dangling, and after: ", sum(collect(values(dangling_inlinks))))
    
    # connect dangling inlinks:

    ap_sec_view = Dict((sec => subset_attachment_potential(ap, sec_view[sec])) for sec in 1:Nnace)

    for (i, ndangling) in dangling_inlinks
        if ndangling > 0
        #     println(i, " ", ndangling)
            new_supps = connect_dangling_bonds(i, nodelist_[i][1], ndangling, Nnace, am2, sec_view, ap_sec_view, supp_view)
            for j in new_supps
                if !((j,i) in edgelist_)
                    addLink!(j,i,nodelist_,edgelist_)
    #             else
    #                 println("Link ($j, $i) already in edgelist_.")
                end
            end
        end
    end
    
    # remove again nodes without links, in case some of the new nodes didn't manage to connect to anyone.
    remove_nodes_wo_links!(nodelist_, Dict{Int64,Int64}(), Dict{Int64,Int64}())
    nothing
end

function run_sim(nl, el, pterm_, prepl_, p_node_exit_, nsp, Nnace, am2, am2_T, fda_, new_nodes_fct_, Nt)
    N_t = zeros(Nt)
    L_t = zeros(Nt)
    for ts in 1:Nt
        println(" --- Timestep $ts --- ")
        timestep!(nl, el, pterm_, prepl_, p_node_exit_, nsp, Nnace, am2, am2_T, fda_, new_nodes_fct_)
        N_t[ts] = length(nl)
        L_t[ts] = length(el)
#         println("Nt = $N_t, Lt = $L_t")
        if (ts%10)==0
            GC.gc()
        end
    end
GC.gc()
return N_t, L_t
end




function nsp_tracker!(new_nodes_sec, removed_nodes_sec,
                      new_links_sec, removed_links_sec,
                      nl_t, el_t, nodelist, edgelist)
    # new_nodes = 0

    for n in keys(nodelist)
        if !(n in keys(nl_t))
    #         new_nodes += 1

            nace_i = nodelist[n][1]
            new_nodes_sec[nace_i] += 1
        end
    end

    # removed_nodes = 0

    for n in keys(nl_t)
        if !(n in keys(nodelist))
    #         removed_nodes += 1

            nace_i = nl_t[n][1]
            removed_nodes_sec[nace_i] += 1
        end
    end 

    # new_links = 0

    for e in edgelist
        if !(e in el_t)
    #         new_links += 1        

            nace_i = nodelist[e[2]][1]  # we're using subject_nace in the other file
            new_links_sec[nace_i] += 1
        end
    end

    # removed_links = 0

    for e in el_t
        if !(e in edgelist)
    #         removed_links += 1

            nace_i = nl_t[e[2]][1]  # we're using subject_nace in the other file
            removed_links_sec[nace_i] += 1
        end
    end
end

# nace_list = unique(collect(keys(nace_consecutive_original_counts)));
function run_sim_w_nsp_tracker(nl, el, pterm_, prepl_, p_node_exit_, nsp, Nnace, am2, am2_T, fda_, new_nodes_fct_, nace_list, Nt)
    N_t = zeros(Nt)
    L_t = zeros(Nt)
    # results = nsp_tracker(nl_t, el_t, nodelist, edgelist)
        
    # make tracking functions:
    new_nodes_sec = Dict([(k,Int(v)) for (k,v) in zip(nace_list, zeros(length(nace_list)))])
    removed_nodes_sec = Dict([(k,Int(v)) for (k,v) in zip(nace_list, zeros(length(nace_list)))])
    new_links_sec = Dict([(k,Int(v)) for (k,v) in zip(nace_list, zeros(length(nace_list)))])
    removed_links_sec = Dict([(k,Int(v)) for (k,v) in zip(nace_list, zeros(length(nace_list)))])

    for ts in 1:Nt
        println(" --- Timestep $ts --- ")
        nl_t = deepcopy(nl)
        el_t = deepcopy(el)
        timestep!(nl, el, pterm_, prepl_, p_node_exit_, nsp, Nnace, am2, am2_T, fda_, new_nodes_fct_)
        N_t[ts] = length(nl)
        L_t[ts] = length(el)
        
        nsp_tracker!(new_nodes_sec, removed_nodes_sec,
                      new_links_sec, removed_links_sec,
                      nl_t, el_t, nl, el)
        
#         println("Nt = $N_t, Lt = $L_t")
        if (ts%10)==0
            GC.gc()
        end
    end
GC.gc()
return N_t, L_t, new_nodes_sec, removed_nodes_sec, new_links_sec, removed_links_sec
end


function save_run(run_name, nodelist, edgelist)
    df_nl = DataFrame([(firm_id=i, k_in=kin, k_out=kout, k_tot=ktot, nace=nace) for (i,(nace, kin, kout, ktot)) in nodelist]);
CSV.write("nodelist_$run_name.csv", df_nl)
    df_el = DataFrame([(supplier_id=i, customer_id=j) for (i,j) in edgelist]);
CSV.write("edgelist_$run_name.csv", df_el)
end

function save_run_ESRI_compatible(run_name, nodelist, edgelist, consecutive_to_nace)
    firm_id_to_consecutive = Dict((firm_id=>i) for (i,firm_id) in enumerate(keys(nodelist)));
    df_nl = DataFrame([(firm_id_consecutive=firm_id_to_consecutive[i]-1, firm_id=i, nace=consecutive_to_nace[sector], k_in=kin, k_out=kout, k_tot=ktot) for (i,(sector, kin, kout, ktot)) in nodelist]);
    sort!(df_nl, [:firm_id_consecutive])
    CSV.write("nodelist_$run_name.csv", df_nl)
    df_el = DataFrame([(supplier=firm_id_to_consecutive[i]-1, buyer=firm_id_to_consecutive[j]-1) for (i,j) in edgelist]);
    CSV.write("edgelist_$run_name.csv", df_el)
end


end


# using Pkg
# Pkg.add("ProgressBars")

# using ProgressBars

# for i in ProgressBar(1:100000) #wrap any iterator
#           #code
#        end