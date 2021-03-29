import LinearAlgebra
import LightGraphs

function to_SimpleGraph_plus_colors(mat::Matrix{Int})
    
    
    @assert all(≥(0),mat) # all labels between 0 and ∞
    @assert LinearAlgebra.issymmetric(mat)
    (m,n) = size(mat)

    g = LightGraphs.SimpleGraph(m)
    e = LightGraphs.edgetype(g)
    edges_color = Dict()
    for i in 1:m
        for j in i+1:n
            if mat[i,j]≠2
                LightGraphs.add_edge!(g,i,j)
                # it seems LightGraphs's algorithm wants both direction for each edge…
                push!(edges_color,e(i,j)=>mat[i,j])
                push!(edges_color,e(j,i)=>mat[i,j])
            end 
        end
    end
    return g, edges_color 

end

function is_isom(mat1,mat2)
    # copied from https://github.com/JuliaGraphs/LightGraphs.jl/blob/ce533729381a6c0b0d8a1eb5f0c3e058833cb55b/src/Structure/isomorphism.jl
    
    g1,c1 = to_SimpleGraph_plus_colors(mat1)
    g2,c2 = to_SimpleGraph_plus_colors(mat2)
    color_rel(e1,e2) = c1[e1] == c2[e2]

    return LightGraphs.Experimental.has_isomorph(g1,g2, LightGraphs.Experimental.VF2(),edge_relation=color_rel)

end

