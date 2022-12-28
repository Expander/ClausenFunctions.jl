@testset "cl1" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Cl1.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test ClausenFunctions.cl1(x) ≈ expected atol=1e-14
        @test ClausenFunctions.cl1(-x) ≈ expected atol=1e-14
        @test ClausenFunctions.cl1(x - 2.0*pi) ≈ expected atol=1e-12
        @test ClausenFunctions.cl1(x + 2.0*pi) ≈ expected atol=1e-12
    end

    @test ClausenFunctions.cl1(1//2) ≈ 0.70358563513784466 atol=1e-14
    @test ClausenFunctions.cl1(1) == ClausenFunctions.cl1(1.0)
    @test ClausenFunctions.cl1(0.0) == Inf
end
