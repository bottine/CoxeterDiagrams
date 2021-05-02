

BS = CoxeterDiagrams.SBitSet

@testset "SBitSet{$N}, to and from" for N in 1:4
    for round in 1:100
        content = Set()
        bs = BS{N}()
        for i in rand(1:N*64,4*N)
            push!(content,i)
            bs = CoxeterDiagrams.push(bs,i)
            @test Set(collect(bs)) == content
        end
        for i in rand(1:N*64,4*N)
            pop!(content,i,nothing)
            bs = CoxeterDiagrams.pop(bs,i)
            @test Set(collect(bs)) == content
        end
    end
end
@testset "SBitSet{$N}, bitwise ops, empty and length" for N in 1:4
    
    data = Vector{Vector{UInt64}}(collect(eachrow(rand(1:64*N,100,4*N))))
    sets = Set.(data)
    bitsets = BS{N}.(data)
    
    @test sets == Set.(collect.(bitsets)) 

    for i in 1:100, j in 1:100
        @test sets[i] ∩ sets[j] == Set(bitsets[i] ∩ bitsets[j])
        @test isempty(sets[i] ∩ sets[j]) == isempty(Set(bitsets[i] ∩ bitsets[j]))
        @test length(sets[i] ∩ sets[j]) == length(Set(bitsets[i] ∩ bitsets[j]))
        @test sets[i] ∪ sets[j] == Set(bitsets[i] ∪ bitsets[j])
        @test isempty(sets[i] ∪ sets[j]) == isempty(Set(bitsets[i] ∪ bitsets[j]))
        @test length(sets[i] ∪ sets[j]) == length(Set(bitsets[i] ∪ bitsets[j]))
        @test setdiff(sets[i],sets[j]) == Set(bitsets[i] ~ bitsets[j])
        @test (sets[i] ⊆ sets[j]) == (bitsets[i] ⊆  bitsets[j])
    end

end
