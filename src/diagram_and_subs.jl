
mutable struct ConnectedInducedSubDiagram
    vertices::SBitSet{4}
    boundary::SBitSet{4}
    type::DiagramType
    degree_sequence::GenDegSeq
    need_to_know_specific_vertices::Bool
    degree_1_vertices::SBitSet{4}
    degree_3::SBitSet{4}
end
CISD = ConnectedInducedSubDiagram

function ConnectedInducedSubDiagram(vertices::SBitSet{4},boundary::SBitSet{4},type::DiagramType)
    ConnectedInducedSubDiagram(vertices,boundary,type,GenDegSeq(GenDeg[]),false,SBitSet{4}(),SBitSet{4}())
end

function ConnectedInducedSubDiagram(v::Int)
    ConnectedInducedSubDiagram(SBitSet{4}(v),SBitSet{4}(v),DT_a,GenDegSeq([empty_deg]),SBitSet{4}(),SBitSet{4}())
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



function connected_diagram_type(VS::SBitSet{4},D::Array{Int,2}; only_sporadic::Bool=false, deg_seq_and_assoc=build_deg_seq_and_associated_data(VS,D))
    
    @assert true "The diagram here is assumed connected. maybe this deserves a check"

   
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
    elseif n≥6 && ds == deg_seq_D(n)  && length(extremities ∩ center_neighbors) ≥ 4
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
    das.D = [das.D; [v;2]']
   
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
    
    _extend!__all_extensions(das,new_vertex,singleton_v,empty_deg,CISD(singleton_v,boundary_v,DT_a,GenDegSeq([empty_deg]),true,SBitSet{4}(),SBitSet{4}()),1,0,new_spherical,new_affine)
   
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
                    
                    if isempty(piece.vertices∩current.vertices) &&  piece.boundary∩current.vertices == singleton_v
                         
                        
                        new_vertices = piece.vertices ∪ current.vertices
                        new_boundary = ((piece.boundary ∩ ~current.vertices) ∪ (current.boundary ∩ ~piece.vertices))

                        new_edges = piece.vertices ∩ current.boundary # all neighbors of v in the new piece we're considering
                        new_edges_col = collect(new_edges)
                        num_new_edges = length(new_edges)
                        new_card = current_card+card
                        new_num_pieces = current_num_pieces + 1
                        

                        if num_new_edges > 2
                            continue
                        elseif num_new_edges == 2 # special case: the only legal way here is if we're closing a diagram of type DT_a to one of type DT_A
                            piece.type == DT_a || continue
                            current_card == 1 || continue
                            (das.D[v,new_edges_col[1]] == 3 && das.D[v,new_edges_col[2]] == 3) || continue
                            new_edges == piece.degree_1_vertices || continue
                            
                            # special case, so we push and are done
                            push!(new_affine[new_card-1],CISD(new_vertices,new_boundary,DT_A))
                        else
                            
                            new_edges_col = collect(new_edges)


                            # only one new edge.
                            u = new_edges_col[1]
                            singleton_u = SBitSet{4}(u)

                            old_deg_seq_u = 0
                            for w in piece.vertices
                                if w ≠ u && das.D[u,w] ≠ 2
                                    old_deg_seq_u = push_label(old_deg_seq_u,das.D[u,w]) 
                                end
                            end

                            l = das.D[v,u]

                            #println("only one new edge")
                            #println("it's $u with label $l")
                            
                            new_deg_seq_v = push_label(deg_seq_v,l) 
                            new_deg_seq_u = push_label(old_deg_seq_u,l)
                            
                            isnothing(new_deg_seq_v) && continue
                            isnothing(new_deg_seq_u) && continue

                            new_gen_deg_seq = current.degree_sequence + piece.degree_sequence
                            update!(new_gen_deg_seq,deg_seq_v,new_deg_seq_v)
                            update!(new_gen_deg_seq,old_deg_seq_u,new_deg_seq_u)
                           
                            new_degree_1 = SBitSet{4}()
                            new_degree_3 = SBitSet{4}()
                            new_degree_1_n = SBitSet{4}()
                            new_degree_3_n = SBitSet{4}()
                            new_need_to_know_specific_vertices = true

                            neighbors(w) = SBitSet{4}([i for i in new_vertices if i≠w && das.D[w,i]== 3])
                           

                            if  (l == 3 || l == 4)  &&
                                current.need_to_know_specific_vertices && 
                                piece.need_to_know_specific_vertices &&
                                ( length(piece.degree_3) + length(current.degree_3) ≤ 2) && 
                                deg_seq_v ∈ (0,1,2,5) && old_deg_seq_u ∈ (0,1,2,5) 

                                new_degree_3 = current.degree_3 ∪ piece.degree_3
                                new_deg_seq_v == 3 && (new_degree_3 = new_degree_3 ∪ singleton_v)
                                new_deg_seq_u == 3 && (new_degree_3 = new_degree_3 ∪ singleton_u)
                                
                                if length(new_degree_3) > 2
                                    
                                    new_need_to_know_specific_vertices = false
                                    
                                else

                                    new_degree_1 = current.degree_1_vertices ∪ piece.degree_1_vertices
                                    l == 3 && old_deg_seq_u == 0 && (new_degree_1 = new_degree_1 ∪ singleton_u)
                                    l == 3 && deg_seq_v == 0 && (new_degree_1 = new_degree_1 ∪ singleton_v)
                                    old_deg_seq_u == 1 && u ∈ new_degree_1 && (new_degree_1 = new_degree_1 ∩ ~ singleton_u)
                                    deg_seq_v == 1 && v ∈ new_degree_1 && (new_degree_1 = new_degree_1 ∩ ~ singleton_v)


                                    new_degree_1_n = SBitSet{4}()
                                    for w in new_degree_1
                                        new_degree_1_n = new_degree_1_n ∪ neighbors(w)
                                    end
                                    new_degree_3_n = SBitSet{4}()
                                    for w in new_degree_3
                                        new_degree_3_n = new_degree_3_n ∪ neighbors(w)
                                    end

                                end 
                                
                            else
                                new_need_to_know_specific_vertices = false
                            end
                             
                            deg_seq_and_assoc = (new_gen_deg_seq,new_degree_1,new_degree_1_n,new_degree_3,new_degree_3_n)
                                 
                            
                            new_type = connected_diagram_type(new_vertices,das.D, only_sporadic=is_sporadic(piece.type) || is_sporadic(current.type), deg_seq_and_assoc=deg_seq_and_assoc)
                            new_deg_1 = new_need_to_know_specific_vertices ? new_degree_1 : SBitSet{4}()
                            new_deg_3 = new_need_to_know_specific_vertices ? new_degree_3 : SBitSet{4}()


                            if !isnothing(new_type)
                                new_cisd = CISD(new_vertices,new_boundary,new_type, new_gen_deg_seq, new_need_to_know_specific_vertices, new_deg_1, new_deg_3)
                                
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

function boundary(das::DiagramAndSubs,vertices::SBitSet{4})
    boundary = SBitSet{4}() 
    neighbors(w) = SBitSet{4}([i for i in 1:size(das.D)[1] if i≠w && das.D[w,i]== 3])
    for w in vertices
        boundary = boundary ∪ neighbors(w)
    end
    boundary = boundary ~ vertices
    return boundary
end

