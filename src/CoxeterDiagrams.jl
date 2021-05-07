
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
    using IterTools
    import Base.show 




    export  DiagramAndSubs, 
            extend!, 
            is_compact, 
            is_finite_volume, 
            is_compact_respectively_finite_volume,  
            is_isom, 
            save_png







    






    include("sbitset2.jl")
    include("diagram_type.jl")
    include("degree_sequence.jl")
    include("isom.jl")
    include("coxiter_io.jl")
    include("util.jl")
    include("diagram_and_subs.jl")
    include("spherical_subdiagram_enumeration.jl")
    include("affine_subdiagram_enumeration.jl")


 
    # Check that all affine subdiagram extend to (at least) a subdiagram of rank das.d-1
    function all_affine_extend_well(das::DiagramAndSubs)
        
        return _all_affine_extend_well__all_extensions(das,SBitSet{4}(),SBitSet{4}(),0)

    end
    function _all_affine_extend_well__all_extensions(
        das::DiagramAndSubs,
        current_vertices::SBitSet{4},
        current_boundary::SBitSet{4},
        current_rank;
        start_rank=1,
        start_idx=1
    )
        
        if current_rank == das.d-1
            return true
        end
        

        has_ext = false 
        @inbounds for piece_rank in start_rank:das.d-1-current_rank
            @inbounds for piece_idx in start_idx:length(das.connected_affine[piece_rank])
                piece = das.connected_affine[piece_rank][piece_idx]
                if  current_rank + piece_rank ≤ das.d-1 &&
                    isempty(piece.vertices∩current_vertices) && 
                    isempty(piece.boundary∩current_vertices) 
                    
                    has_ext = true

                    new_vertices = piece.vertices ∪ current_vertices
                    new_boundary = ((piece.boundary ∩ ~current_vertices) ∪ (current_boundary ∩ ~piece.vertices))
                    new_rank = current_rank + piece_rank 

                    (new_start_rank,new_start_idx) = piece_idx == length(das.connected_affine[piece_rank]) ? (piece_rank+1,1) : (piece_rank,piece_idx+1)
                    if ! _all_affine_extend_well__all_extensions(das,new_vertices,new_boundary,new_rank,start_rank=new_start_rank,start_idx=new_start_idx)
                        return false
                    end
                end
            end  
            start_idx=1
        end

        return has_ext
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
        
        sph_dm = AllSphericalOfRank(das,das.d-1)
        is_empty_sph_dm = true
        for vert in sph_dm
            is_empty_sph_dm = false
            if length(all_spherical_direct_extensions(das,vert)) ≠ 2
                return false
            end
        end
        
        # need to check that there exists a spherical diagram of rank == d.
        # If we are here and there is one of rank d-1, since it has an extension we're good
        # Only remains the case where there is actually no diag of rank d-1
        # more or less degenerate in this case I think
        if is_empty_sph_dm
        #    sph_d = all_spherical_of_rank(das,das.d)
        #    if isempty(sph_d)
                return false
        #    end           
        end

        
        return true
    
    end

    function is_finite_volume(das::DiagramAndSubs; precheck=true)

        # As in Guglielmetti Prop 6.3.1 p. 118 (PhD thesis)
        if precheck && !all_affine_extend_well(das)
            return false
        end
        
        empty_sph_dm = true
        #sph_dm = all_spherical_of_rank(das,das.d-1)
        #for vert in sph_dm 
        for vert in  AllSphericalOfRank(das,das.d-1)
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

        is_empty_sph_dm = true
        sph_dm = AllSphericalOfRank(das,das.d-1)
        for vert in sph_dm
            is_empty_sph_dm = false
            sph_exts = length(all_spherical_direct_extensions(das,vert))
            compact && (compact = (sph_exts == 2))
            if fin_vol 
                aff_exts = length(all_affine_direct_extensions(das,vert))
                fin_vol = (sph_exts + aff_exts == 2)
            end
            !compact && !fin_vol && return (compact,fin_vol)
        end
       
        # more or less degenerate in this case I think
        if is_empty_sph_dm
            return (false,false)
        end
        return (compact,fin_vol)
    end

    function f_vector(das::DiagramAndSubs)
        # TODO: this can be made faster by :
        # * first writing a method `all_spherical` that enumerates all spherical
        # * using this function and adding to the coordinate corresponding to the rank of each diagram we get
        f_vector = Int64[]
        for i in 1:das.d-1 
            push!(f_vector,length(collect(AllSphericalOfRank(das,i))))
        end
        
        push!(f_vector,length(collect(AllSphericalOfRank(das,das.d))))
        f_vector[end] += length(CoxeterDiagrams.all_affine_of_rank(das,das.d-1) |> collect)
        reverse!(f_vector)
        push!(f_vector,1)

        return f_vector
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
        to_neato_io = open(`dot -Tpng`, io, write=true)
        write(to_neato_io,matrix_to_dot(das.D))
        close(to_neato_io)
    end

    function save_png(path, das::DiagramAndSubs)
        open(path,"w") do io
            show(io,MIME"image/png"(),das)
        end
    end

end # end of module




