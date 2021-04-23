using Test
using CSV
using CoxeterDiagrams
using Random
using LinearAlgebra

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

@testset "Compactness/finite volume" begin
    @time begin
    @testset "Known compactness/finite volume values" for row in CSV.Rows("graphs/known_values.csv";comment="#",delim=";",types=[String,Bool,Bool,Float64,String],ignoreemptylines=true)
        @testset "$(row.graph_path)" begin
            println(row.graph_path)
            @test is_compact_respectively_finvol("graphs/"*row.graph_path) == (row.compact,row.finvol)
            #@test is_compact("graphs/"*row.graph_path) == row.compact
        end
    end
    end
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
    end

end
