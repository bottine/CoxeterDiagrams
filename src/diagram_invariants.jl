
     

 

function is_compact(das::DiagramAndSubs)
    
    sph_dm = AllSphericalOfRank(das,das.d-1)
    is_empty_sph_dm = true
    for vert in sph_dm
        is_empty_sph_dm = false
        if number_spherical_direct_extensions_but_at_most_n(das,vert,3) ≠ 2
            return false
        end
    end
    
    # need to check that there exists a spherical diagram of rank == d.
    # If we are here and there is one of rank d-1, since it has an extension we're good
    # Only remains the case where there is actually no diag of rank d-1
    # more or less degenerate in this case I think
    if is_empty_sph_dm
    #    sph_d = all_spherical_of_rank(das,das.d)
    #    if isempty(sph_d)
            return false
    #    end           
    end

    
    return true

end

function is_finite_volume(das::DiagramAndSubs; precheck=false)

    # As in Guglielmetti Prop 6.3.1 p. 118 (PhD thesis)
    if precheck && !all_affine_extend_well(das)
        return false
    end
    
    empty_sph_dm = true
    #sph_dm = all_spherical_of_rank(das,das.d-1)
    #for vert in sph_dm 
    for vert in  AllSphericalOfRank(das,das.d-1)
        empty_sph_dm = false
        sph_exts = number_spherical_direct_extensions_but_at_most_n(das,vert,3)
        aff_exts = number_affine_direct_extensions_but_at_most_n(das,vert,3-sph_exts)
        if !(sph_exts + aff_exts == 2)
            return false
        end
    end
   
    # more or less degenerate in this case I think
    if empty_sph_dm 
        return false
    end
    
    return true
end

function is_compact_respectively_finite_volume(das::DiagramAndSubs)
    
    compact = true
    fin_vol = true

    is_empty_sph_dm = true
    sph_dm = AllSphericalOfRank(das,das.d-1)
    for vert in sph_dm
        is_empty_sph_dm = false
        sph_exts = number_spherical_direct_extensions_but_at_most_n(das,vert,3)
        compact && (compact = (sph_exts == 2))
        if fin_vol 
            aff_exts = number_affine_direct_extensions_but_at_most_n(das,vert,3-sph_exts)
            fin_vol = (sph_exts + aff_exts == 2)
        end
        !compact && !fin_vol && return (compact,fin_vol)
    end
   
    # more or less degenerate in this case I think
    if is_empty_sph_dm
        return (false,false)
    end
    return (compact,fin_vol)
end

function f_vector(das::DiagramAndSubs)
    # * first writing a method `all_spherical` that enumerates all spherical
    # * using this function and adding to the coordinate corresponding to the rank of each diagram we get
    f_vector = [0 for i in 1:das.d]
    for (diagram,rank) in AllSphericalOfRank(das,1,das.d)
        f_vector[rank] += 1
    end
    
    for diagram in AllAffineOfRank(das,das.d-1)
        f_vector[das.d] += 1
    end
    reverse!(f_vector)
    push!(f_vector,1)

    return f_vector
end

