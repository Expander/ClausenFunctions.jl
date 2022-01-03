@testset "cl" begin
    for n in 2:6
        data = open(readdlm, joinpath(@__DIR__, "data", "Cl$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            x        = row[1]
            expected = row[2]
            sgn      = -(-1)^n

            cl  = ClausenFunctions.cl(UInt64(n), x)
            clm = sgn*ClausenFunctions.cl(UInt64(n), -x)

            # println("Cl($(n), $(x)): $(cl) == $(expected)")

            @test cl == clm
            @test is_equal(cl, expected, 1e-10)
        end

    end
end
