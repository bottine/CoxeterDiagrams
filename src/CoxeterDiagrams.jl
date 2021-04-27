
# TODO
#
# * memoize connected_diagram_type by hand
# * make computation of build_deg_seq_and_associated_data more efficient by using previously constructed deg_seqs of the previous components?
# * clean up structure
# * back from bitsets to lists?
#
#
# * be efficient enough to run on the ./graphs/18-vinbxx examples

module CoxeterDiagrams

    
    #import Base.push!, Base.length, Base.copy

    using StaticArrays
    using Memoize
    using Debugger
    using ResumableFunctions
    import Base.show 


    include("sbitset.jl")
    include("diagram_type.jl")
    include("degree_sequence.jl")
    include("isom.jl")
    include("coxiter_io.jl")


    export  DiagramAndSubs, 
            extend!, 
            is_compact, 
            is_finite_volume, 
            is_compact_respectively_finite_volume,  
            is_isom, 
            save_png


    mutable struct ConnectedInducedSubDiagram
        vertices::SBitSet{4}
        boundary::SBitSet{4}
        type::DiagramType
        degree_sequence::GenDegSeq
        need_to_know_specific_vertices::Bool
        degree_1_vertices::Vector{Int}
        degree_1_neighbors::Vector{Int}
        degree_3_neighbors::Vector{Int}
    end
    CISD = ConnectedInducedSubDiagram

    function ConnectedInducedSubDiagram(vertices::SBitSet{4},boundary::SBitSet{4},type::DiagramType)
        ConnectedInducedSubDiagram(vertices,boundary,type,GenDegSeq(GenDeg[]),false,[],[],[])
    end

    function ConnectedInducedSubDiagram(v::Int)
        ConnectedInducedSubDiagram(SBitSet{4}(v),SBitSet{4}(v),DT_a,GenDegSeq([empty_deg]),[],[],[])
    end
    function Base.:(==)(a::CISD,b::CISD)
        # TODO? add the das in CISD so that we can only compare CISDs when they lie in the same DAS?
        a.vertices == b.vertices && a.boundary == b.boundary && a.type == b.type 
    end

    function rank(cisd::ConnectedInducedSubDiagram)
        if is_affine(cisd.type)
            return length(cisd.vertices) - 1
        elseif is_spherical(cisd.type)
            return length(cisd.vertices)
        else
            @assert false "Unreachable: An irreducible diagram is either affine or spherical."
        end
    end

    mutable struct DiagramAndSubs
        D::Array{Int,(2)}                               # = Coxeter matrix
        d::Int                                          # = dimension
        connected_spherical::Vector{Vector{CISD}}               # connected_spherical contains the spherical CISDs 
        connected_affine::Vector{Vector{CISD}}                  # connected_affine contains the affine CISDs
    end

    function is_isom(das1::DiagramAndSubs,das2::DiagramAndSubs)
        return is_isom(das1.D,das2.D) 
    end

    

    function connected_diagram_type(VS::SBitSet{4},D::Array{Int,2}; only_sporadic::Bool=false)
        
        @assert true "The diagram here is assumed connected. maybe this deserves a check"


        deg_seq_and_assoc = build_deg_seq_and_associated_data(VS,D)
       
        if deg_seq_and_assoc === nothing
            return nothing
        end
        

        joined = nothing
        if  length(VS) ≤ 9
            joined = connected_sporadic_diagram_type(VS,D,deg_seq_and_assoc) 
        end
        if joined === nothing && !only_sporadic
            joined = connected_non_sporadic_diagram_type(VS,D,deg_seq_and_assoc)
        end
     


        return joined
    end

    function build_deg_seq_and_associated_data(VS::SBitSet{4},D::Array{Int,2})

        
        @assert true "The diagram here is assumed connected. maybe this deserves a check"

        @debug "build_deg_seq_and_associated_data(…)"
        @debug "VS is $(collect(VS))"

        VSC = length(VS)
       
        deg1 = SBitSet{4}()
        deg1_neigh = SBitSet{4}()
        deg3 = SBitSet{4}()
        deg3_neigh = SBitSet{4}()

        deg_seqs::GenDegSeq = GenDegSeq(Vector{GenDeg}()) 
        for v in VS

            @debug "looking at $v"

            deg_v = empty_deg 
            simple_neighbors_v = SBitSet{4}()
            for u in VS if u ≠ v
                if D[u,v] == 3
                    simple_neighbors_v = simple_neighbors_v | SBitSet{4}(u)
                end
                if D[u,v] ≠ 2
                    deg_v = push_label(deg_v,D[u,v])
                    if deg_v === nothing
                        return nothing
                    end
                end
            end end
            if deg_v == short_vec_to_deg([3,3,3])
                deg3 = deg3 | SBitSet{4}(v)
                deg3_neigh = deg3_neigh | simple_neighbors_v
            elseif deg_v == short_vec_to_deg([3])
                deg1 = deg1 | SBitSet{4}(v)
                deg1_neigh = deg1_neigh | simple_neighbors_v
            end
            # TODO add early exits here
            push!(deg_seqs, deg_v)
            
        end    
        ds = deg_seqs
        

        
        @debug "deg_seq is $ds versus"
        @debug "           $deg_seq_f4"
        @debug "end of build_deg_seq…"


        center::SBitSet{4} = deg3
        center_neighbors::SBitSet{4} = deg3_neigh 
        extremities::SBitSet{4} = deg1 
        extremities_neighbors::SBitSet{4} = deg1_neigh 
        
        return (ds, deg1, deg1_neigh, deg3, deg3_neigh)

    end

    function connected_non_sporadic_diagram_type(VS::SBitSet{4},D::Array{Int,2},deg_seq_and_assoc)
        
        @assert true "The diagram here is assumed connected. maybe this deserves a check"
        
        (ds, extremities, extremities_neighbors, center, center_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 
        
        vertices = VS
        n = length(vertices)

        if false
            @assert false "+++"
        elseif ds == deg_seq_a1
            return DT_a
        elseif n≥2 && ds == deg_seq_a(n)
            return DT_a
        
        elseif ds == deg_seq_b2
            return DT_b
        elseif ds == deg_seq_b3
            return DT_b
        elseif n≥4 && ds == deg_seq_b(n)
            return DT_b
        
        elseif  n≥4 && ds == deg_seq_d(n)  && length(extremities ∩ center_neighbors) ≥ 2
            return DT_d

        elseif n≥3 && ds == deg_seq_A(n)    
            return DT_A

        elseif ds == deg_seq_B4
            return DT_B
        elseif n≥5 && ds == deg_seq_B(n) && length(extremities ∩ center_neighbors) == 2
            return DT_B

        elseif ds == deg_seq_C3
            return DT_C
        elseif n≥4 && ds == deg_seq_C(n)
            return DT_C

        elseif ds == deg_seq_D5 
            return DT_D
        elseif n≥6 && ds == deg_seq_D(n) && length(extremities ∩ center_neighbors) ≥ 4
            return DT_D
        else
            return nothing

        end    
    end

    function connected_sporadic_diagram_type(VS::SBitSet{4},D::Array{Int,2},deg_seq_and_assoc)

        @assert true "The diagram here is assumed connected. maybe this deserves a check"


        (ds, extremities, extremities_neighbors, center, center_neighbors) = deg_seq_and_assoc # build_deg_seq_and_associated_data(VS,D) 
        


        if false
            @assert false "For alignment's sake"
          
        # les sporadiques **pas** de type DT_e ou DT_E
        elseif ds == deg_seq_f4
            return DT_f4
        elseif ds == deg_seq_F4
            return DT_F4
        elseif ds == deg_seq_h2
            return DT_h2
        elseif ds == deg_seq_h3
            return DT_h3
        elseif ds == deg_seq_h4
            return DT_h4
        elseif ds == deg_seq_g2
            return DT_g2
        elseif ds == deg_seq_G2
            return DT_G2
        elseif ds == deg_seq_I∞
            return DT_I∞
        elseif length(ds) == 2  &&
            big_label(ds.content[1]) ≠ nothing &&
            big_label(ds.content[1]) ≥ 7 &&
            ds == deg_seq_i(big_label(ds.content[1]))
            return DT_in


        elseif ds == deg_seq_e6 && length(center_neighbors∩extremities) == 1    
            return DT_e6
        elseif ds == deg_seq_e7 && length(center_neighbors∩extremities) == 1   
            return DT_e7
        elseif ds == deg_seq_e8 && length(center_neighbors∩extremities) == 1 && length(extremities_neighbors ∩ center_neighbors) == 1 
            return DT_e8
        elseif ds == deg_seq_E6 && isempty(center_neighbors∩extremities)    
            return DT_E6
        elseif ds == deg_seq_E7 && length(center_neighbors∩extremities) == 1 && isempty(extremities_neighbors ∩ center_neighbors) 
            return DT_E7
        elseif ds == deg_seq_E8 && length(center_neighbors∩extremities) == 1 && length(extremities_neighbors ∩ center_neighbors) == 1 
            return DT_E8
        end
        
        return nothing

    end



    function extend!(das::DiagramAndSubs, v::Array{Int,1})

        n = length(v)
        @assert size(das.D) == (n,n) "Need $v to have length exactly the number of vertices already present, i.e. $(size(das.D)[1])"
        
        # Extend D with v
        das.D = [das.D v]
        das.D = [das.D; [v;1]']
       
        new_vertex = n+1
        singleton_v = SBitSet{4}(new_vertex)
        boundary_v = SBitSet{4}([i for i in 1:n if das.D[i,new_vertex]≠2])
        new_spherical = [Vector{CISD}() for i in 1:das.d]
        new_affine    = [Vector{CISD}() for i in 1:das.d]
        
        # Add v to the neighbors of preceding connected graphs
        for cisd in Iterators.flatten((Iterators.flatten(das.connected_spherical),Iterators.flatten(das.connected_affine)))
            if !isempty(boundary_v ∩ cisd.vertices)
                cisd.boundary = cisd.boundary ∪ singleton_v
            end
        end
        
        _extend!__all_extensions(das,new_vertex,singleton_v,empty_deg,CISD(singleton_v,boundary_v,DT_a,GenDegSeq([empty_deg]),true,[],[],[]),1,0,new_spherical,new_affine)
       
        for i in 1:das.d
            append!(das.connected_spherical[i],new_spherical[i])
            append!(das.connected_affine[i], new_affine[i])
        end
       
    end
    # compute all possible CISD extensions of `current` and put them in `new_spherical`/`new_affine`
    function _extend!__all_extensions(
        das::DiagramAndSubs,
        v::Int,
        singleton_v,
        deg_seq_v::GenDeg,
        current::CISD,
        current_card::Int,
        current_num_pieces::Int,
        new_spherical,
        new_affine;
        start_card=1,
        start_piece_idx=1
    )
        
        if current_card > das.d
            return
        end
       
        #println("")
        #println("Extending for vertex $v")
        #println("currently $current")
         
            
        if is_spherical(current.type)
            #spherical, so pushed there
            #@assert current_cisd ∉ new_spherical[current_card]
            push!(new_spherical[current_card],current)
            
            if current_num_pieces < 3
                #since spherical, can (probably) be extended:
                @inbounds for card in start_card:length(das.connected_spherical)-current_card+1 # maybe can remove the +1, but it's probably needed for the affine diagrams
                    @inbounds for piece_idx in start_piece_idx:length(das.connected_spherical[card])
                        piece = das.connected_spherical[card][piece_idx]
                        #println("")
                        #println("Extending for vertex $v")
                        #println("currently $current")
                        #println("piec is   $piece")
                        if isempty(piece.vertices∩current.vertices) &&  piece.boundary∩current.vertices == singleton_v
                            
                            #println("intersect nicely OK")
                            
                            new_vertices = piece.vertices ∪ current.vertices
                            new_boundary = ((piece.boundary ∩ ~current.vertices) ∪ (current.boundary ∩ ~piece.vertices))

                            new_edges = collect(piece.vertices ∩ current.boundary) # all neighbors of v in the new piece we're considering
                            num_new_edges = length(new_edges)
                            new_card = current_card+card
                            new_num_pieces = current_num_pieces + 1
                            
                            #println("new edges are $new_edges")
                            #println("num new edges are $num_new_edges")

                            if num_new_edges > 2
                                continue
                            elseif num_new_edges == 2 # special case: the only legal way here is if we're closing a diagram of type DT_a to one of type DT_A
                                piece.type == DT_a || continue
                                current_card == 1 || continue
                                (das.D[v,new_edges[1]] == 3 && das.D[v,new_edges[2]] == 3) || continue
                                new_edges == piece.degree_1_vertices || continue
                                
                                # special case, so we push and are done
                                push!(new_affine[new_card-1],CISD(new_vertices,new_boundary,DT_A))
                            else
                                    
                                # only one new edge.
                                u = new_edges[1]
                                old_deg_seq_u = short_vec_to_deg(das.D[u,filter(≠(u),collect(piece.vertices))])
                                l = das.D[v,u]

                                #println("only one new edge")
                                #println("it's $u with label $l")
                                
                                new_deg_seq_v = push_label(deg_seq_v,l) 
                                new_deg_seq_u = push_label(old_deg_seq_u,l)
                                

                                new_gen_deg_seq = current.degree_sequence + piece.degree_sequence
                                !isnothing(new_deg_seq_v) && update!(new_gen_deg_seq,deg_seq_v,new_deg_seq_v)
                                !isnothing(new_deg_seq_u) && update!(new_gen_deg_seq,old_deg_seq_u,new_deg_seq_u)
                               
                                new_degree_1 = Int[]
                                new_degree_1_n = Int[]
                                new_degree_3_n = Int[]
                                if current.need_to_know_specific_vertices && l == 3 && false 

                                end
                                

                                
                                deg_seq_and_assoc = build_deg_seq_and_associated_data(new_vertices,das.D)
                                if deg_seq_and_assoc ≠ nothing

                                    new_degree_1 = deg_seq_and_assoc[2] |> collect
                                    new_degree_1_n = deg_seq_and_assoc[3] |> collect
                                    new_degree_3_n = deg_seq_and_assoc[5] |> collect

                                    @assert deg_seq_and_assoc[1] == new_gen_deg_seq
                                end

                                new_type = connected_diagram_type(new_vertices,das.D)

                                if !isnothing(new_type)
                                    new_cisd = CISD(new_vertices,new_boundary,new_type, new_gen_deg_seq, current.need_to_know_specific_vertices, new_degree_1, new_degree_1_n, new_degree_3_n)
                                    
                                    (new_start_card,new_start_piece_idx) = piece_idx == length(das.connected_spherical[card]) ? (card+1,1) : (card,piece_idx+1)

                                    #println("continuing with $new_cisd")
                                    _extend!__all_extensions(das,v,singleton_v,new_deg_seq_v,new_cisd,new_card,new_num_pieces,new_spherical,new_affine,start_card=new_start_card,start_piece_idx=new_start_piece_idx)
                                end

                               
                            end

                        end
                    end
                    start_piece_idx = 1
                end               
            end

        elseif is_affine(current.type)
            #affine, so pushed there
            #@assert current_cisd ∉ new_affine[current_card-1]
            push!(new_affine[current_card-1],current)
        else
            @assert false "Diagram either affine or spherical"
        end
            # `current` is invalid (thus cannot be extended to something valid, bye)

        
    end

    


    function DiagramAndSubs(dimension::Int)
        D = reshape(UInt16[],(0,0))
        connected_spherical = [Vector{CISD}() for i in 1:dimension]
        connected_affine    = [Vector{CISD}() for i in 1:dimension]
        return DiagramAndSubs(D,dimension,connected_spherical,connected_affine)
    end

    function DiagramAndSubs(M::Array{Int,2},dimension::Int)
       
        n = size(M,1)
        @assert size(M) == (n,n) "M must be square"
        @assert M == M' "M must be symmetric"
        @assert all(l ≥ 0 for l in M) "M must have non-negative entries"

        das = DiagramAndSubs(dimension)
        for i in 1:n
            extend!(das,M[i,1:i-1])
        end
        return das
    end


    function all_spherical_of_rank(das::DiagramAndSubs,n::Int)
        
        @assert n ≤ das.d
        diagrams_go_here = SBitSet{4}[]#Tuple{SBitSet{4},SBitSet{4}}[]
        _all_spherical_of_rank__all_extensions(das,n,SBitSet{4}(),SBitSet{4}(),0,diagrams_go_here)
        return diagrams_go_here

    end
    function _all_spherical_of_rank__all_extensions(
        das::DiagramAndSubs,
        n::Int,
        current_vertices::SBitSet{4},
        current_boundary::SBitSet{4},
        current_rank::Int,
        diagrams_go_here::Vector{SBitSet{4}};
        start_rank=1,
        start_idx=1
    )
        #@assert current_rank == length(current_vertices)
        
        if current_rank == n
            push!(diagrams_go_here,current_vertices)
            return
        else

            @inbounds for piece_rank in start_rank:n-current_rank
                @inbounds for piece_idx in start_idx:length(das.connected_spherical[piece_rank])
                    piece = das.connected_spherical[piece_rank][piece_idx]
                    #@assert piece_rank == length(piece.vertices)
                    if  current_rank + piece_rank ≤ n &&
                        isempty(piece.vertices ∩ current_vertices) && 
                        isempty(piece.boundary ∩ current_vertices) 
                        
                        new_vertices = piece.vertices ∪ current_vertices
                        new_boundary = ((piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~piece.vertices))
                        new_card = current_rank + piece_rank

                        (new_start_rank,new_start_idx) = piece_idx == length(das.connected_spherical[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)
                        _all_spherical_of_rank__all_extensions(das,n,new_vertices,new_boundary,new_card, diagrams_go_here,start_rank=new_start_rank,start_idx=new_start_idx)
                    end
                end
                start_idx=1
            end 
        end

    end


    function all_affine_of_rank(das::DiagramAndSubs,n::Int)
        
        diagrams_go_here = SBitSet{4}[]#Tuple{SBitSet{4},SBitSet{4}}[]
        _all_affine_of_rank__all_extensions(das,n,SBitSet{4}(),SBitSet{4}(),0,diagrams_go_here)
        return diagrams_go_here

    end
    function _all_affine_of_rank__all_extensions(
        das::DiagramAndSubs,
        n::Int,
        current_vertices::SBitSet{4},
        current_boundary::SBitSet{4},
        current_rank,
        diagrams_go_here::Vector{SBitSet{4}};
        start_rank=1,
        start_idx=1
    )
        
        if current_rank == n
            push!(diagrams_go_here,current_vertices)
            return
        end

        
        @inbounds for piece_rank in start_rank:n-current_rank
            @inbounds for piece_idx in start_idx:length(das.connected_affine[piece_rank])
                piece = das.connected_affine[piece_rank][piece_idx]
                if  current_rank + piece_rank ≤ n &&
                    isempty(piece.vertices∩current_vertices) && 
                    isempty(piece.boundary∩current_vertices) 
                    
                    new_vertices = piece.vertices ∪ current_vertices
                    new_boundary = ((piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~piece.vertices))
                    new_rank = current_rank + piece_rank 

                    (new_start_rank,new_start_idx) = piece_idx == length(das.connected_affine[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)
                    _all_affine_of_rank__all_extensions(das,n,new_vertices,new_boundary,new_rank ,diagrams_go_here,start_rank=new_start_rank,start_idx=new_start_idx)
                end
            end  
            start_idx=1
        end
    end



    
    function all_spherical_direct_extensions(das::DiagramAndSubs,vertices::SBitSet{4})
        
        extensions = SBitSet{4}[]
        for piece in Iterators.flatten(das.connected_spherical)
            if length(piece.vertices ∩ ~vertices) == 1 && isempty(piece.boundary ∩ vertices)
                push!(extensions, piece.vertices ∪ vertices) 
            end
        end
        return extensions
    end

    
    function all_affine_direct_extensions(das::DiagramAndSubs,vertices::SBitSet{4})
        
        diagrams_go_here = SBitSet{4}[]#Tuple{SBitSet{4},SBitSet{4}}[]
        _all_affine_direct_extensions__all_extensions(das,SBitSet{4}(),SBitSet{4}(),vertices,diagrams_go_here)
        return diagrams_go_here
    
    end
    function _all_affine_direct_extensions__all_extensions(
        das::DiagramAndSubs,
        current_vertices::SBitSet{4},
        current_boundary::SBitSet{4},
        remaining_vertices::SBitSet{4},
        diagrams_go_here::Vector{SBitSet{4}};
        start_rank=1,
        start_idx=1
    )
        
        if isempty(remaining_vertices) 
            push!(diagrams_go_here,current_vertices)
            return
        end


        @inbounds for piece_rank in start_rank:length(das.connected_affine)
            @inbounds for piece_idx in start_idx:length(das.connected_affine[piece_rank])
                piece = das.connected_affine[piece_rank][piece_idx]
                if  isempty(piece.vertices∩current_vertices) && 
                    isempty(piece.boundary∩current_vertices) &&
                    length(piece.vertices ∩ ~remaining_vertices) == 1 &&
                    isempty(piece.boundary ∩ remaining_vertices)
                    
                    new_vertices = piece.vertices ∪ current_vertices
                    new_boundary = ((piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~piece.vertices))
                    new_remaining_vertices = remaining_vertices ∩ ~piece.vertices
                    
                    (new_start_rank,new_start_idx) = piece_idx == length(das.connected_affine[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)
                    _all_affine_direct_extensions__all_extensions(das,new_vertices,new_boundary,new_remaining_vertices,diagrams_go_here,start_rank=new_start_rank,start_idx=new_start_idx)
                end
            end
            start_idx = 1
        end               
    end
   

    function is_compact(das::DiagramAndSubs)
        
        sph_dm = all_spherical_of_rank(das,das.d-1)
        for vert in sph_dm
            if length(all_spherical_direct_extensions(das,vert)) ≠ 2
                return false
            end
        end
        
        # need to check that there exists a spherical diagram of rank == d.
        # If we are here and there is one of rank d-1, since it has an extension we're good
        # Only remains the case where there is actually no diag of rank d-1
        # more or less degenerate in this case I think
        if isempty(sph_dm)
        #    sph_d = all_spherical_of_rank(das,das.d)
        #    if isempty(sph_d)
                return false
        #    end           
        end

        
        return true
    
    end

    function is_finite_volume(das::DiagramAndSubs)
        
        empty_sph_dm = true
        #sph_dm = all_spherical_of_rank(das,das.d-1)
        #for vert in sph_dm 
        for vert in  all_spherical_of_rank(das,das.d-1)
            empty_sph_dm = false
            sph_exts = length(all_spherical_direct_extensions(das,vert))
            aff_exts = length(all_affine_direct_extensions(das,vert))
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

        sph_dm = all_spherical_of_rank(das,das.d-1)
        for vert in sph_dm
            sph_exts = length(all_spherical_direct_extensions(das,vert))
            compact && (compact = (sph_exts == 2))
            if fin_vol 
                aff_exts = length(all_affine_direct_extensions(das,vert))
                fin_vol = (sph_exts + aff_exts == 2)
            end
            !compact && !fin_vol && return (compact,fin_vol)
        end
       
        # more or less degenerate in this case I think
        if isempty(sph_dm)
            return (false,false)
        end
        return (compact,fin_vol)
    end

    function is_compact_respectively_finite_volume(path::String)
        
        s = open(path) do file
            read(file, String)
        end

        ret = gug_coxiter_to_matrix(s)
        if ret === nothing
            println("Failed reading $path")
        else
            (D, rank) = ret
            if D === nothing || rank === nothing
                println("Error reading file probably")
            else
                das = DiagramAndSubs(D,rank)
                (compact,fin_vol) = is_compact_respectively_finite_volume(das)
                return (compact,fin_vol) 
            end
        end


    end


   #= 

    function check_all_graphs(sub_directory="")
        

        for (root, dirs, files) in walkdir("../graphs/"*sub_directory)
            for path in joinpath.(root, files)
                if endswith(path,".coxiter")
                    println("path: $path")
                    println(is_compact_respectively_finvol(path))
                end
            end
        end

    end


    function check_some_graphs()

        @time begin
            for path in ["../graphs/13-mcl11.coxiter"] 
                    println("path: $path")
                    println(is_compact_respectively_finvol(path))
            end
            for (root, dirs, files) in walkdir("../graphs/simplices")
                for path in joinpath.(root, files)
                    if endswith(path,".coxiter")
                        println("path: $path")
                        println(is_compact_respectively_finvol(path))
                    end
                end
            end
        end

    end

    =#

    function matrix_to_dot(Mat)
        
        (m,n) = size(Mat)
        
        @assert m == n

        dot_string = """
    strict graph {
        layout=neato
        """

        for i in 1:m
            dot_string = dot_string * "\tnode [shape=circle,label=\"$i\"] $i;\n"
        end


        function label_to_edge_type(k)
            if k == 1
                #return "[style=dotted,label=$k,weight=0]"
                return "[style=dotted,weight=0]"
            elseif k == 0
                #return "[penwidth=3,label=$k,weight=2]"
                return "[penwidth=3,weight=2]"
            elseif k == 3
                #return "[label=$k]"
                return ""
            elseif k > 3
                #return "[color = \"" * "black:invis:"^(k-3) * ":black\",label=$k]"
                return "[color = \"" * "black:invis:"^(k-3) * ":black\"]"
            end
        end

        for i in 1:m
            for j in i+1:n
                if Mat[i,j] ≠ 2
                    dot_string = dot_string * "\t$i -- $j" * label_to_edge_type(Mat[i,j]) * ";\n"
                end
            end
        end

        dot_string = dot_string * "}"


        return dot_string

    end

    

    function Base.show(io, ::MIME"image/png", das::DiagramAndSubs)
        to_neato_io = open(`neato -Tpng`, io, write=true)
        write(to_neato_io,matrix_to_dot(das.D))
        close(to_neato_io)
    end

    function save_png(path, das::DiagramAndSubs)
        open(path,"w") do io
            show(io,MIME"image/png"(),das)
        end
    end

end # end of module




