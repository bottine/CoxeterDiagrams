

mutable struct Stack{T,N}
    stack::MVector{N,T}
    idx::Int
end

Stack{T,N}() where {T,N} = Stack{T,N}(MVector{N,T}(undef),0)

function Base.push!(s::Stack{T,N},t::T) where {T,N}
    @assert s.idx < N
    s.idx += 1
    s.stack[s.idx] = t
    return s
end

function Stack{T,N}(t::T) where {T,N}
    s = Stack{T,N}(MVector{N,T}(undef),0)
    push!(s,t)
    return s
end

function Base.pop!(s::Stack{T,N}) where {T,N}
    @assert s.idx > 0
    s.idx -= 1
    return s.stack[s.idx+1]
end

function Base.isempty(s::Stack{T,N}) where {T,N}
    s.idx == 0
end

