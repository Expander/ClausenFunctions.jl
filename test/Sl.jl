@testset "sl" begin
    Sl1(x) = 0.5*(pi - x)
    Sl2(x) = pi^2/6 - pi/2*x + x^2/4
    Sl3(x) = pi^2/6*x - pi/4*x^2 + x^3/12
    Sl4(x) = pi^4/90 - pi^2/12*x^2 + pi/12*x^3 - x^4/48

    n = 10
    x_min = 0.0
    x_max = 2.0*pi
    data = (x_max - x_min)*rand(Float64, n) + x_min*ones(n)

    for x in data
        @test ≈(ClausenFunctions.sl(1, x), Sl1(x), atol=1e-14)
        @test ≈(ClausenFunctions.sl(2, x), Sl2(x), atol=1e-14)
        @test ≈(ClausenFunctions.sl(3, x), Sl3(x), atol=1e-14)
        @test ≈(ClausenFunctions.sl(4, x), Sl4(x), atol=1e-13)
    end

    for n in vcat(collect(1:31), [1000, 1001])
        data = open(readdlm, joinpath(@__DIR__, "data", "Sl$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            x        = row[1]
            expected = row[2]

            @test ClausenFunctions.sl(n, x) == (-1)^n*ClausenFunctions.sl(n, -x)
            @test ClausenFunctions.sl(n, x) ≈ expected atol=1e-14
        end
    end

    @test_throws DomainError ClausenFunctions.sl(0, 1.0)
    @test_throws DomainError ClausenFunctions.sl(-1, 1.0)
    @test_throws DomainError ClausenFunctions.sl(-2, 1.0)
end
