
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

function AllAffineOfRank(das,min,max)
    return AllAffineOfRank(das,min,max,SBitSet{4}())
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
                
                @tassert piece_rank == length(piece.vertices) - 1
                @tassert current_rank + piece_rank ≤ max
                
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



function all_affine_extend_well(das)
    affine = SBitSet{4}[]
    affine_rank_dm = SBitSet{4}[]

    for (diag,rank) in AllAffineOfRank(das,1,das.d-1)
        push!(affine,diag)
        rank == das.d-1 && push!(affine_rank_dm, diag) 
    end
    
    for diag in affine

        if !any(diag ⊆ diag_dm for diag_dm in affine_rank_dm)
            return false
        end
    end
    return true
    
end


function all_affine_extend_well_bad(das)
   
    stack = Stack{Tuple{SBitSet{4},SBitSet{4},Int,Int,Int,Bool},20}((SBitSet{4}(),SBitSet{4}(),0,1,1,false)) # Because we assume das.d ≤ 20 always!
    
    @label killbillvol2
    while !isempty(stack)
        
        println("-------------------------------------")
        println("stack is:")
        for s in stack.stack[1:stack.idx]
            println("$(s[1]) $(s[6])")
        end
        println("-------------------------------------")
        println("looking at")

        (current_vertices,current_boundary,current_rank,start_rank,start_idx,needs_extension_yet) = pop!(stack)
        println("$current_vertices $needs_extension_yet ($start_rank, $start_idx)")
        println("------------------------------------")
        extends = !needs_extension_yet 

        @inbounds for piece_rank in start_rank:das.d-1-current_rank
            @inbounds for piece_idx in start_idx:length(das.connected_affine[piece_rank])
                piece = das.connected_affine[piece_rank][piece_idx]
                
                println("can the following extend?")
                println(piece.vertices)

                @tassert piece_rank == length(piece.vertices) - 1
                
                if  isempty(piece.vertices ∩ current_vertices) && 
                    isempty(piece.boundary ∩ current_vertices) 
                   
                    println("extendedable with:")
                    println("$(piece.vertices)")
                    extends = true

                    new_vertices = piece.vertices ∪ current_vertices
                    new_boundary = ((piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~piece.vertices))
                    new_rank = current_rank + piece_rank

                    (new_start_rank,new_start_idx) = piece_idx == length(das.connected_affine[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)

                    if new_rank == das.d-1

                        push!(stack, (current_vertices,current_boundary,current_rank,new_start_rank,new_start_idx,false))
                    else
                        push!(stack, (current_vertices,current_boundary,current_rank,new_start_rank,new_start_idx,needs_extension_yet))
                        push!(stack, (new_vertices,new_boundary,new_rank,1,1,true))
                    end
                    @goto killbillvol2
                end
            end
            start_idx=1
        end
        
        !extends && return false

    end
    
    return true

end





struct AllAffineDirectExtensions
    das::DiagramAndSubs
    base::SBitSet{4}
end

function Base.iterate(a::AllAffineDirectExtensions)
    state = Stack{Tuple{SBitSet{4},SBitSet{4},SBitSet{4},Int,Int},20}((SBitSet{4}(),SBitSet{4}(),a.base,1,1)) # Because we assume das.d ≤ 20 always!
    
    return iterate(a,state)
end

function Base.iterate(a::AllAffineDirectExtensions,state)
   
    das = a.das
    stack = state

    @label killbillvol2
    while !isempty(stack)
    
        (current_vertices,current_boundary,remaining_vertices,start_rank,start_idx) = pop!(stack)
        if isempty(remaining_vertices) 
            return (current_vertices,stack)
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
                    push!(stack,(current_vertices,current_boundary,remaining_vertices,new_start_rank,new_start_idx))
                    push!(stack,(new_vertices,new_boundary,new_remaining_vertices,new_start_rank,new_start_idx))
                    @goto killbillvol2
                end
            end
            start_idx = 1
        end  
    end
end

Base.IteratorEltype(::Type{AllAffineDirectExtensions}) = Base.HasEltype() 
Base.eltype(::Type{AllAffineDirectExtensions}) = SBitSet{4}
Base.IteratorSize(::Type{AllAffineDirectExtensions}) = Base.SizeUnknown()

function all_affine_direct_extensions(das::DiagramAndSubs,vertices::SBitSet{4})
    
    return AllAffineDirectExtensions(das,vertices)

end

function number_affine_direct_extensions_but_at_most_n(das::DiagramAndSubs,vertices::SBitSet{4},n)
    num_exts = 0
    for ext in all_affine_direct_extensions(das,vertices)
        num_exts += 1
        num_exts ≥ n && return num_exts
    end
    return num_exts
end

#=

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
=#
