@testset "cl" begin
    for n in [2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16, 1000, 1001, 1_000_000]
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
            @test cl ≈ expected atol=1e-14
        end

    end

    @test_throws DomainError ClausenFunctions.cl(1, 1.0)
    @test_throws DomainError ClausenFunctions.cl(0, 1.0)
    @test_throws DomainError ClausenFunctions.cl(-1, 1.0)
end
