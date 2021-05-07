

@testset "Finite volume precheck negative implies not finite volume, and not finite volume implies not compact." begin
    @testset "Known diagrams" for row in CSV.Rows("../graphs/known_values.csv";
                                                                          comment="#",
                                                                          delim=";",
                                                                          types=[String,Int,Int,Int,Int,Int,Bool,Bool,Float64,String],ignoreemptylines=true)
        @testset "$(row.graph_path)" begin

            path = row.graph_path 
            s = open("../graphs/"*path) do file
                read(file, String)
            end

            ret = CoxeterDiagrams.gug_coxiter_to_matrix(s)
            if ret === nothing
                println("Failed reading $path")
            else
                (D, rank) = ret
                if D === nothing || rank === nothing
                    println("Error reading file probably")
                else

                    println("$path ($(size(D)[1]))") 
                    das = DiagramAndSubs(rank)
                    for i in 1:size(D)[1]
                        extend!(das,D[i,1:i-1])
                        println("$i")
                        precheck = CoxeterDiagrams.all_affine_extend_well(das)
                        finvol = is_finite_volume(das,precheck=false)
                        compact = is_compact(das)

                        #@test precheck3 == precheck2 
                        @test precheck ≥ finvol
                        @test finvol ≥ compact 
                    end

                    println()

                end
            end


        end
    end
end

