using Test
using CSV
using CoxeterDiagrams
using Random
using LinearAlgebra

CD = CoxeterDiagrams
CISD = CD.ConnectedInducedSubDiagram
SBS = CD.SBitSet{4}


include("unit/sbitset2.jl")

@testset "Degree Sequences" begin
    
end

@testset "Some randomized Coxeter matrix isomorphism checks" begin
   
    max_card = 20
    max_coeff = 5
    @testset "Cardinality $card" for card in 3:max_card
        
        for run in 10*(max_card-card)
           
            M = LinearAlgebra.tril!(rand(0:max_coeff,card,card)) # Take a random matrix, make it triangular
            M = M + M'                                           # and symmetric
            
            σ = shuffle(1:card)                                  # Find a random permutation
            σM = M[σ,σ]                                          # permute the matrix accordingly
            @test CoxeterDiagrams.is_isom(M,σM)
            
            if card>2
                # change one edge
                mut_i = rand(1:card-1)
                mut_j = rand(mut_i+1:card)
                mut_v = rand(1:max_coeff)
                mut_M = copy(M)
                mut_M[mut_i,mut_j] = M[mut_i,mut_j]+mut_v
                mut_M[mut_j,mut_i] = M[mut_j,mut_i]+mut_v

                @test !CoxeterDiagrams.is_isom(M,mut_M)
                @test !CoxeterDiagrams.is_isom(σM,mut_M)
            end

            N = LinearAlgebra.tril!(rand(0:max_coeff,card,card))
            N = N + N'

            # different degree sequences
            M_sorted = sort([sort(c) for c in eachcol(M)])
            N_sorted = sort([sort(c) for c in eachcol(N)])
            if M_sorted ≠ N_sorted
                @test !CoxeterDiagrams.is_isom(M,N)
            end

        end
    end

end


@testset "One specific example" begin
    mat = [2 3 3 2 2 2 2;
           3 2 3 2 2 2 2;
           3 3 2 4 2 2 2;
           2 2 4 2 3 3 2;
           2 2 2 3 2 3 2;
           2 2 2 3 3 2 6;
           2 2 2 2 2 6 2]
    das = DiagramAndSubs(mat,7)
    
    affine_should_be = CISD[
        CISD(SBS([1, 2, 3]), SBS([4]), CD.DT_A), 
        CISD(SBS([4, 5, 6]), SBS([3, 7]), CD.DT_A), 
        CISD(SBS([4, 6, 7]), SBS([3, 5]), CD.DT_G2), 
        CISD(SBS([5, 6, 7]), SBS([4]), CD.DT_G2)
    ]

    @test length(Set(vcat(das.connected_affine...))) == length(vcat(das.connected_affine...))
    @test length(affine_should_be) == length(vcat(das.connected_affine...)) 

    spherical_should_be = CISD[
         CISD(SBS([1]), SBS([2, 3]), CD.DT_a),
         CISD(SBS([2]), SBS([1, 3]), CD.DT_a),
         CISD(SBS([1, 2]), SBS([3]), CD.DT_a),
         CISD(SBS([3]), SBS([1, 2, 4]), CD.DT_a),
         CISD(SBS([1, 3]), SBS([2, 4]), CD.DT_a),
         CISD(SBS([2, 3]), SBS([1, 4]), CD.DT_a),
         CISD(SBS([4]), SBS([3, 5, 6]), CD.DT_a),
         CISD(SBS([3, 4]), SBS([1, 2, 5, 6]), CD.DT_b),
         CISD(SBS([1, 3, 4]), SBS([2, 5, 6]), CD.DT_b),
         CISD(SBS([2, 3, 4]), SBS([1, 5, 6]), CD.DT_b),
         CISD(SBS([5]), SBS([4, 6]), CD.DT_a),
         CISD(SBS([4, 5]), SBS([3, 6]), CD.DT_a),
         CISD(SBS([3, 4, 5]), SBS([1, 2, 6]), CD.DT_b),
         CISD(SBS([1, 3, 4, 5]), SBS([2, 6]), CD.DT_f4),
         CISD(SBS([2, 3, 4, 5]), SBS([1, 6]), CD.DT_f4),
         CISD(SBS([6]), SBS([4, 5, 7]), CD.DT_a),
         CISD(SBS([4, 6]), SBS([3, 5, 7]), CD.DT_a),
         CISD(SBS([3, 4, 6]), SBS([1, 2, 5, 7]), CD.DT_b),
         CISD(SBS([1, 3, 4, 6]), SBS([2, 5, 7]), CD.DT_f4),
         CISD(SBS([2, 3, 4, 6]), SBS([1, 5, 7]), CD.DT_f4),
         CISD(SBS([5, 6]), SBS([4, 7]), CD.DT_a),
         CISD(SBS([7]), SBS([6]), CD.DT_a),
         CISD(SBS([6, 7]), SBS([4, 5]), CD.DT_g2)
    ]
  
    @test length(Set(vcat(das.connected_spherical...))) == length(vcat(das.connected_spherical...))
    @test length(spherical_should_be) == length(vcat(das.connected_spherical...))
   
    @test [length(CD.all_affine_of_rank(das,i) |> collect) for i in 1:5] == [0, 4, 0, 1, 0]
   
    @test [length(collect(CoxeterDiagrams.all_spherical_of_rank(das,i))) for i in 1:6] == [7, 21, 31, 21, 3, 0] 

end





@testset "More random matrices" begin
    for rank in 2:13, size in 5:20, round in 1:10
        M = rand([0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,4,5,6,7,8],size,size); M = tril(M); M = M + M'
        das = DiagramAndSubs(M,rank)
        compact,fin_vol = is_compact_respectively_finite_volume(das)

        @testset "is_compact_respectively_finite_volume agrees with is_compact and is_finite_volume" begin
            @test compact == is_compact(das)
            @test fin_vol == is_finite_volume(das)
        end

        
        @testset "Computing direct extensions is the same as filtering over the candidate extensions" for i in 1:rank-1 
            all_sph_i = CoxeterDiagrams.all_spherical_of_rank(das,i) |> collect
            all_sph_ip = CoxeterDiagrams.all_spherical_of_rank(das,i+1) |> collect
            all_aff_i = CoxeterDiagrams.all_affine_of_rank(das,i) |> collect
            for sph in all_sph_i
                ext_sph_1 = CoxeterDiagrams.all_spherical_direct_extensions(das,sph)
                ext_sph_2 = [sph2 for sph2 in all_sph_ip if sph ⊆ sph2]
                @test Set(ext_sph_1) == Set(ext_sph_2) 
                ext_aff_1 = CoxeterDiagrams.all_affine_direct_extensions(das,sph) 
                ext_aff_2 = [aff2 for aff2 in all_aff_i if sph ⊆ aff2]
                @test Set(ext_aff_1) == Set(ext_aff_2) 
            end
        end

    end
end

@testset "f-vector must ±sum to 1" begin
    @time begin
    @testset "On known diagrams" for row in CSV.Rows("../graphs/known_values.csv";
                                                                          comment="#",
                                                                          delim=";",
                                                                          missingstring="",
                                                                          types=[String,Int,Int,Int,Int,Int,Bool,Bool,Float64,String],ignoreemptylines=true)
        if !occursin("drop",row.graph_path) 
        @testset "$(row.graph_path)" begin
            println("$(row.graph_path)")
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
                    das = DiagramAndSubs(D,rank)
                    euler = 0
                    f_vector = CoxeterDiagrams.f_vector(das)
                    for (ip,f_i) in enumerate(vcat([1],f_vector))
                        i = ip-1
                        euler += (-1)^i * f_i
                    end
                    @test euler == 0
                end
            end
        end
        end
    end
    end
end


@testset "Compactness/finite volume" begin
    @time begin
    @testset "Known compactness/finite volume values" for row in CSV.Rows("../graphs/known_values.csv";
                                                                          comment="#",
                                                                          delim=";",
                                                                          missingstring="",
                                                                          types=[String,Int,Int,Int,Int,Int,Bool,Bool,Float64,String],ignoreemptylines=true)
        @testset "$(row.graph_path)" begin
            println(row.graph_path)
            @test is_compact_respectively_finite_volume("../graphs/"*row.graph_path) == (row.compact,row.finvol)
        end
    end
    end
end


@testset "Right number of irreducible spherical/affine subdiagrams" begin
    @testset "Known compactness/finite volume values" for row in CSV.Rows("../graphs/known_values.csv";
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
                    das = DiagramAndSubs(D,rank)
                    if !ismissing(row.num_irred_sph) 
                        @test sum(das.connected_spherical .|> length) == row.num_irred_sph
                    end
                    if !ismissing(row.num_irred_aff) 
                        @test  sum(das.connected_affine .|> length) == row.num_irred_aff
                    end
                    if !ismissing(row.num_sph_rank_dm) 
                        @test  CoxeterDiagrams.all_spherical_of_rank(das,das.d-1) |> collect |> length  == row.num_sph_rank_dm
                    end
                    if !ismissing(row.num_sph_rank_d) 
                        @test  CoxeterDiagrams.all_spherical_of_rank(das,das.d)  |> collect |> length == row.num_sph_rank_d
                    end 
                    if !ismissing(row.num_aff_rank_dm) 
                        @test  CoxeterDiagrams.all_affine_of_rank(das,das.d-1)  |> collect |> length == row.num_aff_rank_dm
                    end
                    #@test CoxeterDiagrams.is_compact_respectively_finite_volume(das) == (is_compact(das),is_finite_volume(das))
                    #println("> ", length(all_spherical_of_rank(das,das.d-1)), ", ", length(all_spherical_of_rank(das,das.d)), ", ", length(all_affine_of_rank(das,das.d-1)))
                    #dump_das(das;range=nothing)
                    #compact = is_compact(das)
                    #fin_vol = is_finite_volume(das)
                    #@assert (compact,fin_vol) == (is_compact(das),is_finite_volume(das))
                end
            end


        end
    end
end

