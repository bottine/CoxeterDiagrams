
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
    using Debugger
    import Base.show 


    include("sbitset.jl")
    include("diagram_type.jl")
    include("degree_sequence.jl")
    include("isom.jl")
    include("coxiter_io.jl")


    export build_diagram_and_subs, extend!, is_compact, is_finite_volume, is_compact_finite_volume, is_fin_vol, is_compact_respectively_finvol, is_isom, save_png


    mutable struct ConnectedInducedSubDiagram
        vertices::SBitSet{4}
        boundary::SBitSet{4}
        type::DiagramType
        #degree_sequence::GenDegSeq
        #degree_1_vertices::Vector{UInt16}
        #degree_3_vertex::Vector{UInt16}
    end
    CISD = ConnectedInducedSubDiagram

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
        connected_spherical::Vector{CISD}               # connected_spherical contains the spherical CISDs 
        connected_affine::Vector{CISD}                  # connected_affine contains the affine CISDs
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
        new_spherical = Vector{CISD}()
        new_affine = Vector{CISD}()


        
        # Add v to the neighbors of preceding connected graphs
        for cisd in Iterators.flatten((das.connected_spherical,das.connected_affine))
            if !isempty(boundary_v ∩ cisd.vertices)
                cisd.boundary = cisd.boundary ∪ singleton_v
            end
        end


        # compute all possible CISD extensions of `current` and put them in `new_spherical`/`new_affine`
        function all_extensions(current_vertices::SBitSet{4},current_boundary::SBitSet{4};look_after=1)
            
            if length(current_vertices) > das.d + 1
                return
            end

            # Is `current` a valid connected diagram?
            current_type = connected_diagram_type(current_vertices,das.D)
            if !isnothing(current_type)

                # if here, it is a valid connected diagram, which we "wrap" in a nice type
                current_cisd = CISD(current_vertices,current_boundary,current_type)
                
                if is_spherical(current_type)
                    #spherical, so pushed there
                    @assert current_cisd ∉ new_spherical; push!(new_spherical,current_cisd)
                    
                    #since spherical, can (probably) be extended:
                    for (idx,new_piece) in enumerate(das.connected_spherical[look_after:end])
                        if isempty(new_piece.vertices∩current_vertices) &&  new_piece.boundary∩current_vertices == singleton_v
                            new_vertices = new_piece.vertices ∪ current_vertices
                            new_boundary = ((new_piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~new_piece.vertices))
                            all_extensions(new_vertices,new_boundary,look_after=look_after+idx)
                        end
                    end               
                elseif is_affine(current_type)
                    #affine, so pushed there
                    @assert current_cisd ∉ new_affine; push!(new_affine,current_cisd)
                else
                    @assert false "Diagram either affine or spherical"
                end
            else
                # `current` is invalid (thus cannot be extended to something valid, bye)
            end
        end
        
        all_extensions(singleton_v,boundary_v)
        
        append!(das.connected_spherical,new_spherical)
        append!(das.connected_affine, new_affine)
       
    end


    function DiagramAndSubs(dimension::Int)
        D = reshape(UInt16[],(0,0))
        connected_spherical = Vector{CISD}()
        connected_affine    = Vector{CISD}()
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
        
        diagrams_go_here = SBitSet{4}[]#Tuple{SBitSet{4},SBitSet{4}}[]
        function all_extensions(current_vertices::SBitSet{4},current_boundary::SBitSet{4};look_after=1)
            
            if length(current_vertices) == n
                push!(diagrams_go_here,current_vertices)
                #push!(diagrams_go_here,(current_vertices,current_boundary))
                return
            end

            for (idx,new_piece) in enumerate(das.connected_spherical[look_after:end])
                if  length(new_piece.vertices) + length(current_vertices) ≤ n &&
                    isempty(new_piece.vertices∩current_vertices) && 
                    isempty(new_piece.boundary∩current_vertices) 
                    
                    new_vertices = new_piece.vertices ∪ current_vertices
                    new_boundary = ((new_piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~new_piece.vertices))
                    all_extensions(new_vertices,new_boundary,look_after=look_after+idx)
                end
            end               
        end
         
        all_extensions(SBitSet{4}(),SBitSet{4}())

        return diagrams_go_here
    end

    function all_affine_of_rank(das::DiagramAndSubs,n::Int)
        
        diagrams_go_here = SBitSet{4}[]#Tuple{SBitSet{4},SBitSet{4}}[]
        function all_extensions(current_vertices::SBitSet{4},current_boundary::SBitSet{4},current_rank;look_after=1)
            
            if current_rank == n
                push!(diagrams_go_here,current_vertices)
                #push!(diagrams_go_here,(current_vertices,current_boundary))
                return
            end

            for (idx,new_piece) in enumerate(das.connected_affine[look_after:end])
                if  current_rank + length(new_piece.vertices) - 1 ≤ n &&
                    isempty(new_piece.vertices∩current_vertices) && 
                    isempty(new_piece.boundary∩current_vertices) 
                    
                    new_vertices = new_piece.vertices ∪ current_vertices
                    new_boundary = ((new_piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~new_piece.vertices))
                    all_extensions(new_vertices,new_boundary,current_rank + length(new_piece.vertices) - 1,look_after=look_after+idx)
                end
            end               
        end
         
        all_extensions(SBitSet{4}(),SBitSet{4}(),0)

        return diagrams_go_here
    end




    
    function all_spherical_direct_extensions(das::DiagramAndSubs,vertices::SBitSet{4})
        
        extensions = SBitSet{4}[]
        for piece in das.connected_spherical
            if length(piece.vertices ∩ ~vertices) == 1 && isempty(piece.boundary ∩ vertices)
                push!(extensions, piece.vertices ∪ vertices) 
            end
        end
        return extensions
    end

    
    function all_affine_direct_extensions(das::DiagramAndSubs,vertices::SBitSet{4})
        

        diagrams_go_here = SBitSet{4}[]#Tuple{SBitSet{4},SBitSet{4}}[]
        function all_extensions(
            current_vertices::SBitSet{4},
            current_boundary::SBitSet{4},
            remaining_vertices::SBitSet{4};
            look_after=1
        )
            
            if isempty(remaining_vertices) 
                push!(diagrams_go_here,current_vertices)
                return
            end

            for (idx,new_piece) in enumerate(das.connected_affine[look_after:end])
                if  isempty(new_piece.vertices∩current_vertices) && 
                    isempty(new_piece.boundary∩current_vertices) &&
                    length(new_piece.vertices ∩ ~remaining_vertices) == 1 &&
                    isempty(new_piece.boundary ∩ remaining_vertices)
                    
                    new_vertices = new_piece.vertices ∪ current_vertices
                    new_boundary = ((new_piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~new_piece.vertices))
                    new_remaining_vertices = remaining_vertices ∩ ~new_piece.vertices
                    all_extensions(new_vertices,new_boundary,new_remaining_vertices,look_after=look_after+idx)
                end
            end               
        end
         
        all_extensions(SBitSet{4}(),SBitSet{4}(),vertices)


        return diagrams_go_here
    end
    

    function is_compact(das::DiagramAndSubs)
        length(all_spherical_of_rank(das,das.d)) > 0 && 
        all(
            length(all_spherical_direct_extensions(das,vert)) == 2 for 
            vert in all_spherical_of_rank(das,das.d-1)
        )
    end

    function is_finvol(das::DiagramAndSubs)

        sph_d = all_spherical_of_rank(das,das.d)
        sph_dm = all_spherical_of_rank(das,das.d-1)
        aff_dm = all_affine_of_rank(das,das.d-1)

        if length(all_spherical_of_rank(das,das.d)) == 0 && length(all_affine_of_rank(das,das.d-1)) == 0
            return false
        end

        for vert in all_spherical_of_rank(das,das.d-1)
            sph_exts = length(all_spherical_direct_extensions(das,vert))
            aff_exts = length(all_affine_direct_extensions(das,vert))
            #aff_exts = length(filter(aff -> vert ⊆ aff, aff_dm))
            if !(sph_exts + aff_exts == 2)
                return false
            end
        end
        return true
    end

    # ##########################
    #
    # Misc user-facing functions
    #
    # ##########################
    #

    function is_fin_vol(D::Array{Int,2},dimension::Int)
        
        das = build_diagram_and_subs(D,dimension)
        return is_finite_volume(das)

    end

    function is_compact(path::String)
        
        # Get file content in s as in https://en.wikibooks.org/wiki/Introducing_Julia/Working_with_text_files
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
                return is_compact(das) 
            end
        end
    end
     
    function is_compact_respectively_finvol(path::String)
        
        # Get file content in s as in https://en.wikibooks.org/wiki/Introducing_Julia/Working_with_text_files
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
                #dump_das(das;range=nothing)
                compact = is_compact(das)
                fin_vol = is_finvol(das)
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
        node [shape=point];
        """

        function label_to_edge_type(k)
            if k == 1
                return "[style=dotted,label=$k,weight=0]"
            elseif k == 0
                return "[penwidth=3,label=$k,weight=2]"
            elseif k == 3
                return "[label=$k]"
            elseif k > 3
                return "[color = \"" * "black:invis:"^(k-3) * ":black\",label=$k]"
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




