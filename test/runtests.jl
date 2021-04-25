using Test
using CSV
using CoxeterDiagrams
using Random
using LinearAlgebra

CD = CoxeterDiagrams
CISD = CD.ConnectedInducedSubDiagram
SBS = CD.SBitSet{4}

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
   
    @test [length(CD.all_affine_of_rank(das,i)) for i in 1:5] == [0, 4, 0, 1, 0]
    #@test length(CoxeterDiagrams.all_affine_of_rank(das,1)) == 0
    #@test length(CoxeterDiagrams.all_affine_of_rank(das,2)) == 4
    #@test length(CoxeterDiagrams.all_affine_of_rank(das,3)) == 0
    #@test length(CoxeterDiagrams.all_affine_of_rank(das,4)) == 1
    #@test length(CoxeterDiagrams.all_affine_of_rank(das,5)) == 0
   
    @test [length(CoxeterDiagrams.all_spherical_of_rank(das,i)) for i in 1:6] == [7, 21, 31, 21, 3, 0] 

end




@testset "Compactness/finite volume for some randomly generated matrices (checks that the code agrees with old versions of itself)" begin
    
    list = Tuple{Int64, Matrix{Int64}, Bool, Bool}[
         (3, [8 1 2 0 2 6; 1 2 0 0 1 0; 2 0 4 2 0 2; 0 0 2 8 2 6; 2 1 0 2 4 2; 6 0 2 6 2 12], 0, 1),
		 (2, [10 1 0 2 6; 1 0 2 6 1; 0 2 0 0 2; 2 6 0 4 0; 6 1 2 0 0], 1, 0),
		 (4, [4 1 3 2 2 2; 1 4 1 1 2 1; 3 1 4 3 3 2; 2 1 3 4 2 3; 2 2 3 2 2 2; 2 1 2 3 2 4], 1, 1),
		 (3, [2 0 0 0 1; 0 4 2 2 2; 0 2 0 2 2; 0 2 2 4 5; 1 2 2 5 4], 1, 1),
		 (3, [2 2 2 1 2 2; 2 4 1 2 2 2; 2 1 2 3 5 2; 1 2 3 0 2 2; 2 2 5 2 0 1; 2 2 2 2 1 0], 1, 1),
		 (3, [4 2 0 2 2; 2 0 2 0 2; 0 2 2 2 0; 2 0 2 4 2; 2 2 0 2 10], 0, 1),
		 (3, [4 0 3 2 2; 0 4 2 2 2; 3 2 8 3 2; 2 2 3 6 0; 2 2 2 0 2], 0, 1),
		 (3, [2 6 1 0 2 2; 6 4 2 1 2 2; 1 2 4 2 3 6; 0 1 2 2 2 2; 2 2 3 2 4 0; 2 2 6 2 0 6], 0, 1),
		 (2, [8 0 2 0 2; 0 4 1 2 2; 2 1 4 2 1; 0 2 2 6 1; 2 2 1 1 0], 1, 0),
		 (3, [0 2 1 2 3; 2 2 1 2 2; 1 1 2 0 1; 2 2 0 0 2; 3 2 1 2 8], 1, 1),
		 (3, [0 2 3 2 0; 2 6 2 0 2; 3 2 8 3 2; 2 0 3 4 2; 0 2 2 2 0], 0, 1),
		 (2, [4 1 0 2 2; 1 4 2 0 2; 0 2 6 2 1; 2 0 2 4 1; 2 2 1 1 2], 1, 0),
		 (3, [4 6 2 2 2; 6 4 2 3 2; 2 2 4 2 0; 2 3 2 4 2; 2 2 0 2 4], 1, 0),
		 (4, [4 5 1 2 3 2; 5 2 0 2 2 2; 1 0 0 4 1 2; 2 2 4 6 2 4; 3 2 1 2 2 3; 2 2 2 4 3 10], 1, 1),
		 (2, [4 2 6 1 1; 2 10 1 0 1; 6 1 4 1 4; 1 0 1 4 0; 1 1 4 0 8], 0, 1),
		 (2, [2 0 2 2 0; 0 4 1 5 2; 2 1 2 1 2; 2 5 1 4 0; 0 2 2 0 10], 1, 0),
		 (4, [2 2 2 2 2; 2 0 2 2 2; 2 2 0 2 2; 2 2 2 0 2; 2 2 2 2 12], 1, 1),
         (4, [4 3 2 3 2 1; 3 0 2 2 2 1; 2 2 0 2 2 1; 3 2 2 0 2 1; 2 2 2 2 4 2; 1 1 1 1 2 2], 1, 1), 
         (4, [2 2 2 2 1 2 2; 2 4 2 2 6 2 0; 2 2 4 2 0 4 0; 2 2 2 10 4 2 0; 1 6 0 4 2 1 6; 2 2 4 2 1 6 1; 2 0 0 0 6 1 0], 1, 1), 
         (4, [8 2 2 4 2; 2 4 2 2 2; 2 2 2 2 2; 4 2 2 12 2; 2 2 2 2 12], 1, 1), 
         (5, [2 2 3 2 2 4 3; 2 4 2 2 3 4 2; 3 2 0 2 2 0 2; 2 2 2 2 2 0 3; 2 3 2 2 4 0 2; 4 4 0 0 0 2 2; 3 2 2 3 2 2 4], 1, 1)
    ]  

    for (rank,mat,compact,fin_vol) in list
        das = CoxeterDiagrams.DiagramAndSubs(mat,rank)
        @test compact == CoxeterDiagrams.is_compact(das)
        @test fin_vol == CoxeterDiagrams.is_finite_volume(das)
        @test (compact,fin_vol) == CoxeterDiagrams.is_compact_respectively_finite_volume(das)
    end

end

#=
@testset "Right number of irreducible spherical/affine subdiagrams" begin
    @testset "Known compactness/finite volume values" for row in CSV.Rows("../graphs/known_values.csv";
                                                                          comment="#",
                                                                          delim=";",
                                                                          missingstring="",
                                                                          types=[String,Int,Int,Int,Int,Int,Bool,Bool,Float64,String],ignoreemptylines=true)
        @testset "$(row.graph_path)" begin
            println(row.graph_path)
            println(row)

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
                    println(row.num_irred_sph)
                    println(row.num_irred_aff)
                    if !ismissing(row.num_irred_sph) 
                        @test sum(das.connected_spherical .|> length) == row.num_irred_sph
                    end
                    if !ismissing(row.num_irred_aff) 
                        @test  sum(das.connected_affine .|> length) == row.num_irred_aff
                    end
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
=#

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
