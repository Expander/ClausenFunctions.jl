@testset "cl2" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Cl2.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test ClausenFunctions.cl2(x) ≈ expected rtol=1e-14 atol=1e-14
        @test ClausenFunctions.cl2(-x) ≈ -expected rtol=1e-14 atol=1e-14
        @test ClausenFunctions.cl2(x - 2pi) ≈ expected rtol=1e-13 atol=1e-13
        @test ClausenFunctions.cl2(x + 2pi) ≈ expected rtol=1e-13 atol=1e-13

        @test ClausenFunctions.cl2(Float16(x)) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
        @test ClausenFunctions.cl2(Float32(x)) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
    end

    @test ClausenFunctions.cl2(0.0) == 0.0
    @test ClausenFunctions.cl2(pi) == 0.0
    @test ClausenFunctions.cl2(pi/2) ≈ 0.915965594177219015054603514932384110 rtol=1e-14
    @test ClausenFunctions.cl2(1//2) ≈ 0.84831187770367927 rtol=1e-14

    # test handling of negative zero
    @test !signbit(ClausenFunctions.cl2(0.0))
    @test signbit(ClausenFunctions.cl2(-0.0))
end
