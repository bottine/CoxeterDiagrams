

@testset "all_affine_extend_well on random matrices" begin
    for rank in 2:13, size in rank+1:20, round in 1:10
        M = rand([0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,4,5,6,7,8],size,size); M = tril(M); M = M + M'
        das = DiagramAndSubs(M,rank)

        @test CoxeterDiagrams.all_affine_extend_well_safe(das) == CoxeterDiagrams.all_affine_extend_well(das)

    end
end




@testset "all_affine_extend_well on known examples" begin
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
                        @test CoxeterDiagrams.all_affine_extend_well_safe(das) == CoxeterDiagrams.all_affine_extend_well(das)
                    end

                end
            end


        end
    end
end

