#
# Freely inspired from there: https://raw.githubusercontent.com/chethega/StaticArrays.jl/fb0350012f01db4021d60906357e949333ec5d93/src/SBitSet.jl
#



struct SBitSet{N}
    pieces::NTuple{N,UInt64}
end


@inline function SBitSet(m::MBitSet{N}) where N
    SBitSet{N}(m.pieces)
end
@inline function MBitSet(s::SBitSet{N}) where N
    SBitSet{N}(s.pieces)
end


@inline function Base.length(s::SBitSet{N}) where N
    sum(count_ones.(s.pieces))
end

@inline function Base.isempty(s::SBitSet{N}) where N
    all(==(0), s.pieces)
end

@inline function Base.:(==)(s1::SBitSet{N},s2::SBitSet{N}) where N
    s1.pieces == s2.pieces
end

@inline function Base.:(&)(s1::SBitSet{N},s2::SBitSet{N}) where N
    SBitSet{N}((&).(s1.pieces,s2.pieces))
end
Base.:(∩)(s1::SBitSet{N},s2::SBitSet{N}) where N = Base.:(&)(s1,s2)

@inline function Base.:(|)(s1::SBitSet{N},s2::SBitSet{N}) where N
    SBitSet{N}((|).(s1.pieces,s2.pieces))
end
Base.:(∪)(s1::SBitSet{N},s2::SBitSet{N}) where N = Base.:(|)(s1,s2)

@inline function Base.:(~)(s1::SBitSet{N}) where N
    SBitSet{N}((~).(s1.pieces))
end
Base.:(~)(s1::SBitSet{N},s2::SBitSet{N}) where N = s1 ∩ ~s2

@inline function Base.:(⊆)(s1::SBitSet{N},s2::SBitSet{N}) where N
   s1 ∩ s2 == s1 
end

@inline function SBitSet{N}() where N
    SBitSet{N}(ntuple(x->UInt64(0),N))
end

@inline function SBitSet{N}(n::Integer) where N
    @assert 1 ≤ n && n ≤ N*64
    (d,r) = divrem(n,64)
    if r == 0
        (d,r) = (d,64)
    else
        (d,r) = (d+1,r)
    end
    SBitSet{N}(ntuple(x-> x == d ? (UInt64(1) << (r-1)) : UInt64(0),N))
end

@inline function SBitSet{N}(v::Vector) where N
    s = SBitSet{N}()
    for i in v
        s = push(s,i)
    end
    return s
end

@inline function push(s::SBitSet{N},n::Integer) where N
    @assert 1 ≤ n && n ≤ N*64
    (d,r) = divrem(n,64)
    if r == 0
        (d,r) = (d,64)
    else
        (d,r) = (d+1,r)
    end
    SBitSet{N}(ntuple(x-> x == d ? s.pieces[x] | (UInt64(1) << (r-1)) : s.pieces[x],N))
end

@inline function pop(s::SBitSet{N},n::Integer) where N
    @assert 1 ≤ n && n ≤ N*64
    (d,r) = divrem(n,64)
    if r == 0
        (d,r) = (d,64)
    else
        (d,r) = (d+1,r)
    end
    SBitSet{N}(ntuple(x-> x == d ? s.pieces[x] & ~(UInt64(1) << (r-1)) : s.pieces[x],N))
end

@inline function Base.iterate(s::SBitSet{N}) where N
    Base.iterate(s,(UInt64(1),UInt64(1),UInt64(1),s.pieces[1]))
end

@inline function Base.iterate(s::SBitSet{N},(d,r,current,current_block)) where N
    d>N && return nothing
    #@assert (d-1)*64 + r == current
    while d ≤ N
        while current_block ≠ UInt64(0)

            #= semi random to speed up going over holes =#
            while current_block & UInt64(0xFFFF) == 0
                current_block >>= 16
                r += 16
                current += 16
            end 
            while current_block & UInt64(0xF) == 0
                current_block >>= 4
                r += 4
                current += 4
            end
            
            if current_block & UInt64(1) == UInt64(1)
                return (current, (d,r+1,current+1,current_block>>1))
            end
            current+=1
            r+=1
            current_block >>= 1
        end
        d += 1
        r = 1
        current = 64*(d-1)+1
        d > N && return nothing
        current_block = s.pieces[d]
    end
    
    return nothing 
end


Base.eltype(::Type{SBitSet}) = UInt64
