@testset "cl4" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Cl4.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test is_equal(ClausenFunctions.cl4(x), expected, 1e-14)
        @test is_equal(ClausenFunctions.cl4(-x), -expected, 1e-14)
        @test is_equal(ClausenFunctions.cl4(x - 2.0*pi), expected, 1e-13)
        @test is_equal(ClausenFunctions.cl4(x + 2.0*pi), expected, 1e-13)
    end
end
