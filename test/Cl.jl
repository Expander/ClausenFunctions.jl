@testset "cl" begin
    for n in 2:6
        data = open(readdlm, joinpath(@__DIR__, "data", "Cl$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            x        = row[1]
            expected = row[2]
            sgn      = -(-1)^n

            cl  = ClausenFunctions.cl(n, x)
            clm = sgn*ClausenFunctions.cl(n, -x)

            # println("Cl($(n), $(x)): $(cl) == $(expected)")

            @test cl == clm
            @test is_equal(cl, expected, 1e-14)
        end

    end

    @test_throws DomainError ClausenFunctions.cl(1, 1.0)
    @test_throws DomainError ClausenFunctions.cl(0, 1.0)
    @test_throws DomainError ClausenFunctions.cl(-1, 1.0)
end
