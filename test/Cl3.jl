@testset "cl3" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Cl3.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test ClausenFunctions.cl3(x) ≈ expected atol=1e-14
        @test ClausenFunctions.cl3(-x) ≈ expected atol=1e-14
        @test ClausenFunctions.cl3(x - 2.0*pi) ≈ expected atol=1e-13
        @test ClausenFunctions.cl3(x + 2.0*pi) ≈ expected atol=1e-13

        @test ClausenFunctions.cl3(Float16(x)) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
        @test ClausenFunctions.cl3(Float32(x)) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
    end

    @test ClausenFunctions.cl3(1//2) ≈ 0.92769631047023043 atol=1e-14
end
