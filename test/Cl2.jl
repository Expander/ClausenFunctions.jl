@testset "cl2" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Cl2.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test ClausenFunctions.cl2(x) ≈ expected atol=1e-14
        @test ClausenFunctions.cl2(-x) ≈ -expected atol=1e-14
        @test ClausenFunctions.cl2(x - 2.0*pi) ≈ expected atol=1e-13
        @test ClausenFunctions.cl2(x + 2.0*pi) ≈ expected atol=1e-13
    end
end
