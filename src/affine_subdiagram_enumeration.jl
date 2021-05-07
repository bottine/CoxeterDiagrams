
struct AllAffineOfRank
    das::DiagramAndSubs
    min::Int
    max::Int
    base::SBitSet{4}
    function AllAffineOfRank(das::DiagramAndSubs,min::Int,max::Int,base::SBitSet{4})
        @assert 1 ≤ min ≤ max ≤ das.d
        return new(das,min,max,base)
    end
end

function AllAffineOfRank(das,n)
    return IterTools.imap(x->x[1],AllAffineOfRank(das,n,n,SBitSet{4}()))
end

function Base.iterate(a::AllAffineOfRank)
    state = Stack{Tuple{SBitSet{4},SBitSet{4},Int,Int,Int},20}((a.base,boundary(a.das,a.base),length(a.base),1,1)) # Because we assume das.d ≤ 20 always!
    # TODO: If base is a valid diagram, also emit it
    
    return iterate(a,state)
end

function Base.iterate(a::AllAffineOfRank,state)
   
    das = a.das
    min = a.min
    max = a.max
    stack = state

    @label killbillvol2
    while !isempty(stack)

        (current_vertices,current_boundary,current_rank,start_rank,start_idx) = pop!(stack)

        @inbounds for piece_rank in start_rank:max-current_rank
            @inbounds for piece_idx in start_idx:length(das.connected_affine[piece_rank])
                piece = das.connected_affine[piece_rank][piece_idx]
                
                @assert piece_rank == length(piece.vertices) - 1
                @assert current_rank + piece_rank ≤ max
                
                if  isempty(piece.vertices ∩ current_vertices) && 
                    isempty(piece.boundary ∩ current_vertices) 
                    
                    new_vertices = piece.vertices ∪ current_vertices
                    new_boundary = ((piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~piece.vertices))
                    new_rank = current_rank + piece_rank

                    (new_start_rank,new_start_idx) = piece_idx == length(das.connected_affine[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)

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

Base.IteratorEltype(::Type{AllAffineOfRank}) = Base.HasEltype() 
Base.eltype(::Type{AllAffineOfRank}) = Tuple{SBitSet{4},Int}
Base.IteratorSize(::Type{AllAffineOfRank}) = Base.SizeUnknown()

all_affine_of_rank(das::DiagramAndSubs,n::Int) = AllAffineOfRank(das,n) 
all_affine_of_rank(das::DiagramAndSubs,min::Int,max::Int) = AllAffineOfRank(das,min,max,SBitSet{4}())

