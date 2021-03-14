using Test
using CSV
using CoxeterDiagrams

@testset "Compactness/finite volume" begin

    @testset "Known compactness/finite volume values" for row in CSV.Rows("../graphs/known_values.csv";comment="#",delim=";",types=[String,Bool,Bool,Float64,String],ignoreemptylines=true)
        @test is_compact_respectively_finvol("../graphs/"*row.graph_path) == (row.compact,row.finvol)
    end

end
