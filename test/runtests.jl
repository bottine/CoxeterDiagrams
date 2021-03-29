using Test
using CSV
using CoxeterDiagrams
using Random
using LinearAlgebra

@testset "Coxeter Matrix Isomorphism Check" begin
   
    max_card = 20
    max_coeff = 5
    @testset "Cardinality $card" for card in 2:max_card
        
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

    @testset "Known compactness/finite volume values" for row in CSV.Rows("../graphs/known_values.csv";comment="#",delim=";",types=[String,Bool,Bool,Float64,String],ignoreemptylines=true)
        @test is_compact_respectively_finvol("../graphs/"*row.graph_path) == (row.compact,row.finvol)
    end

end
