
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


    # Copied from the ToggleableAsserts source
    module Options

        const options_lock = ReentrantLock()

        assert_enabled() = false

        function toggle_asserts(b::Bool)
            lock(options_lock) do
                @eval Options assert_enabled() = $b
            end
        end

        macro tassert(cond, text=nothing)
            if text==nothing
                assert_stmt = esc(:(@assert $cond))
            else
                assert_stmt = esc(:(@assert $cond $text))
            end
            :(assert_enabled() ? $assert_stmt  : nothing)
        end
                

        export @tassert, toggle_asserts 
    end

    
    using CoxeterDiagrams.Options
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
    include("diagram_invariants.jl")    


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




