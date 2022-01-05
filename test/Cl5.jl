@testset "cl5" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Cl5.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test ClausenFunctions.cl5(x) ≈ expected atol=1e-14
        @test ClausenFunctions.cl5(-x) ≈ expected atol=1e-14
        @test ClausenFunctions.cl5(x - 2.0*pi) ≈ expected atol=1e-13
        @test ClausenFunctions.cl5(x + 2.0*pi) ≈ expected atol=1e-13
    end
end
