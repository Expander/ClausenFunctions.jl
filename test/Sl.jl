@testset "sl" begin
    Sl1(x) = 0.5*(pi - x)
    Sl2(x) = pi^2/6 - pi/2*x + x^2/4
    Sl3(x) = pi^2/6*x - pi/4*x^2 + x^3/12
    Sl4(x) = pi^4/90 - pi^2/12*x^2 + pi/12*x^3 - x^4/48

    n = 10
    x_min = 0.0
    x_max = 2pi
    data = (x_max - x_min)*rand(Float64, n) + x_min*ones(n)

    for x in data
        @test ClausenFunctions.sl(1, x) ≈ Sl1(x) rtol=1e-14 atol=1e-14
        @test ClausenFunctions.sl(2, x) ≈ Sl2(x) rtol=1e-14 atol=1e-14
        @test ClausenFunctions.sl(3, x) ≈ Sl3(x) rtol=1e-14 atol=1e-14
        @test ClausenFunctions.sl(4, x) ≈ Sl4(x) rtol=1e-14 atol=1e-13
    end

    for n in vcat(collect(1:31), [1000, 1001])
        data = open(readdlm, joinpath(@__DIR__, "data", "Sl$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            x        = row[1]
            expected = row[2]

            @test ClausenFunctions.sl(n, x) == (-1)^n*ClausenFunctions.sl(n, -x)
            @test ClausenFunctions.sl(n, x) ≈ expected rtol=1e-14 atol=1e-14
            @test ClausenFunctions.sl(n, Float16(x)) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
            @test ClausenFunctions.sl(n, Float32(x)) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
            @test ClausenFunctions.sl(n, Complex(x)) ≈ expected rtol=1e-14 atol=1e-14
            @test ClausenFunctions.sl(n, Complex(Float16(x))) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
            @test ClausenFunctions.sl(n, Complex(Float32(x))) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
        end
    end

    @test ClausenFunctions.sl(1, 1//2) ≈ 1.3207963267948966 rtol=1e-14
    @test ClausenFunctions.sl(2, 1//2) ≈ 0.92203590345077813 rtol=1e-14
    @test ClausenFunctions.sl(3, 1//2) ≈ 0.63653415924141781 rtol=1e-14
    @test ClausenFunctions.sl(4, 1//2) ≈ 0.90812931549667023 rtol=1e-14
    @test ClausenFunctions.sl(5, 1//2) ≈ 0.51085256423059275 rtol=1e-14
    @test ClausenFunctions.sl(6, 1//2) ≈ 0.88593812938731573 rtol=1e-14

    @test ClausenFunctions.sl(1, Complex(1//2)) ≈ 1.3207963267948966 rtol=1e-14
    @test ClausenFunctions.sl(2, Complex(1//2)) ≈ 0.92203590345077813 rtol=1e-14
    @test ClausenFunctions.sl(3, Complex(1//2)) ≈ 0.63653415924141781 rtol=1e-14
    @test ClausenFunctions.sl(4, Complex(1//2)) ≈ 0.90812931549667023 rtol=1e-14
    @test ClausenFunctions.sl(5, Complex(1//2)) ≈ 0.51085256423059275 rtol=1e-14
    @test ClausenFunctions.sl(6, Complex(1//2)) ≈ 0.88593812938731573 rtol=1e-14

    # @test ClausenFunctions.sl(1, 0.0 + 1.0im) ≈ -0.0413248546129181 - 1.5707963267948966im rtol=1e-14
    # @test ClausenFunctions.sl(1, 0.0 - 1.0im) ≈ -0.0413248546129181 - 1.5707963267948966im rtol=1e-14
    # @test ClausenFunctions.sl(2, 0.0 + 1.0im) ≈ 1.5707963267948966 + 0.9861797794993302im rtol=1e-14
    # @test ClausenFunctions.sl(2, 0.0 - 1.0im) ≈ -1.5707963267948966 - 0.9861797794993302im rtol=1e-14
    # @test ClausenFunctions.sl(2, 1.0 + 1.0im) ≈ 1.4107754938116412 - 0.1044778629291566im rtol=1e-14
    # @test ClausenFunctions.sl(2, 1.0 - 1.0im) ≈ 1.4107754938116412 + 0.1044778629291566im rtol=1e-14

    @test_throws DomainError ClausenFunctions.sl(0, 1.0)
    @test_throws DomainError ClausenFunctions.sl(-1, 1.0)
    @test_throws DomainError ClausenFunctions.sl(-2, 1.0)
end
