import LinearAlgebra
using LightGraphs

function to_SimpleGraph_plus_colors(mat::Matrix{Int})
    
    
    @assert all(≥(0),mat) # all labels between 0 and ∞
    @assert LinearAlgebra.issymmetric(mat)
    (m,n) = size(mat)

    vertices = Any[i for i in 1:n]
    vertices_c = Int[-1 for i in 1:n]
    adj = Any[]
    for i in 1:m
        for j in i+1:n
            if mat[i,j]≠2
                edge = (i,j)
                push!(vertices,edge)
                push!(vertices_c,mat[i,j])
                push!(adj,(i,edge))
                push!(adj,(j,edge))
            end 
        end
    end
    
    g = SimpleGraph(length(vertices))

    for (u,v) in adj
        idx_u = indexin(u,vertices)[1]
        idx_v = indexin(v,vertices)[1]
        add_edge!(g,idx_u,idx_v)
    end

    return g, vertices_c 

end

function is_isom(mat1,mat2)
    # copied from https://github.com/JuliaGraphs/LightGraphs.jl/blob/ce533729381a6c0b0d8a1eb5f0c3e058833cb55b/src/Structure/isomorphism.jl
    
    g1,c1 = to_SimpleGraph_plus_colors(mat1)
    g2,c2 = to_SimpleGraph_plus_colors(mat2)
    color_rel(v1,v2) = c1[v1] == c2[v2]

    return LightGraphs.Experimental.has_isomorph(g1,g2, LightGraphs.Experimental.VF2(),vertex_relation=color_rel)

end

