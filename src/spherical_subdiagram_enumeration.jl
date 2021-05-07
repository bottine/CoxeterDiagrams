


struct AllSphericalOfRank
    das::DiagramAndSubs
    min::Int
    max::Int
    base::SBitSet{4}
    function AllSphericalOfRank(das::DiagramAndSubs,min::Int,max::Int,base::SBitSet{4})
        @assert 1 ≤ min ≤ max ≤ das.d
        return new(das,min,max,base)
    end
end

function AllSphericalOfRank(das,n)
    return IterTools.imap(x->x[1],AllSphericalOfRank(das,n,n,SBitSet{4}()))
end
function AllSphericalOfRank(das,min,max)
    return AllSphericalOfRank(das,min,max,SBitSet{4}())
end

function Base.iterate(a::AllSphericalOfRank)
    state = Stack{Tuple{SBitSet{4},SBitSet{4},Int,Int,Int},20}((a.base,boundary(a.das,a.base),length(a.base),1,1)) # Because we assume das.d ≤ 20 always!
    # TODO: If base is a valid diagram, also emit it
    
    return iterate(a,state)
end

function Base.iterate(a::AllSphericalOfRank,state)
   
    das = a.das
    min = a.min
    max = a.max
    stack = state

    @label killbillvol2
    while !isempty(stack)

        (current_vertices,current_boundary,current_rank,start_rank,start_idx) = pop!(stack)

        @inbounds for piece_rank in start_rank:max-current_rank
            @inbounds for piece_idx in start_idx:length(das.connected_spherical[piece_rank])
                piece = das.connected_spherical[piece_rank][piece_idx]
                

                @tassert piece_rank == length(piece.vertices)
                @tassert current_rank + piece_rank ≤ max
                
                if  isempty(piece.vertices ∩ current_vertices) && 
                    isempty(piece.boundary ∩ current_vertices) 
                    
                    new_vertices = piece.vertices ∪ current_vertices
                    new_boundary = ((piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~piece.vertices))
                    new_rank = current_rank + piece_rank

                    (new_start_rank,new_start_idx) = piece_idx == length(das.connected_spherical[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)

                    push!(stack, (current_vertices,current_boundary,current_rank,new_start_rank,new_start_idx))
                    new_rank ≠ max && push!(stack, (new_vertices,new_boundary,new_rank,new_start_rank,new_start_idx))
                    if min ≤ new_rank
                        return ((new_vertices,new_rank),stack)
                    else
                        @goto killbillvol2
                    end
                end
            end
            start_idx=1
        end 

    end


end

Base.IteratorEltype(::Type{AllSphericalOfRank}) = Base.HasEltype() 
Base.eltype(::Type{AllSphericalOfRank}) = Tuple{SBitSet{4},Int}
Base.IteratorSize(::Type{AllSphericalOfRank}) = Base.SizeUnknown()

all_spherical_of_rank(das::DiagramAndSubs,n::Int) = AllSphericalOfRank(das,n) 
all_spherical_of_rank(das::DiagramAndSubs,min::Int,max::Int) = AllSphericalOfRank(das,min,max,SBitSet{4}())



struct AllSphericalDirectExtensions
    das::DiagramAndSubs
    base::SBitSet{4}
end

function Base.iterate(a::AllSphericalDirectExtensions,state=(1,1))
    start_rank,start_idx = state
    das = a.das
    vertices = a.base

    @inbounds for piece_rank in start_rank:das.d
        @inbounds for piece_idx in start_idx:length(das.connected_spherical[piece_rank])
            piece = das.connected_spherical[piece_rank][piece_idx]
            if length(piece.vertices ∩ ~vertices) == 1 && isempty(piece.boundary ∩ vertices)
                (new_start_rank,new_start_idx) = piece_idx == length(das.connected_spherical[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)
                return (piece.vertices ∪ vertices,(new_start_rank,new_start_idx)) 
            end
        end
        start_idx=1
    end
end

Base.IteratorEltype(::Type{AllSphericalDirectExtensions}) = Base.HasEltype() 
Base.eltype(::Type{AllSphericalDirectExtensions}) = SBitSet{4}
Base.IteratorSize(::Type{AllSphericalDirectExtensions}) = Base.SizeUnknown()

function all_spherical_direct_extensions2(das::DiagramAndSubs,vertices::SBitSet{4})
    return AllSphericalDirectExtensions(das,vertices)
end

function 


function all_spherical_direct_extensions(das::DiagramAndSubs,vertices::SBitSet{4})
    
    extensions = SBitSet{4}[]
    for piece in Iterators.flatten(das.connected_spherical)
        if length(piece.vertices ∩ ~vertices) == 1 && isempty(piece.boundary ∩ vertices)
            push!(extensions, piece.vertices ∪ vertices) 
        end
    end

    @assert extensions == all_spherical_direct_extensions2(das,vertices) |> collect "$extensions vs $(all_spherical_direct_extensions2(das,vertices) |> collect) for $vertices in $(das.D)"

    return extensions
end

