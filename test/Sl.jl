@testset "sl" for TN in (Int8, Int16, Int32, Int64, Int128)
    Sl1(x) = 0.5*(pi - x)
    Sl2(x) = pi^2/6 - pi/2*x + x^2/4
    Sl3(x) = pi^2/6*x - pi/4*x^2 + x^3/12
    Sl4(x) = pi^4/90 - pi^2/12*x^2 + pi/12*x^3 - x^4/48

    n = 10
    x_min = 0.0
    x_max = 2pi
    data = (x_max - x_min)*rand(Float64, n) + x_min*ones(n)

    for x in data
        @test ClausenFunctions.sl(TN(1), x) ≈ Sl1(x) rtol=1e-14 atol=1e-14
        @test ClausenFunctions.sl(TN(2), x) ≈ Sl2(x) rtol=1e-14 atol=1e-14
        @test ClausenFunctions.sl(TN(3), x) ≈ Sl3(x) rtol=1e-14 atol=1e-14
        @test ClausenFunctions.sl(TN(4), x) ≈ Sl4(x) rtol=1e-14 atol=1e-13
    end

    for n in vcat(collect(-2:31), [1000, 1001])
        (n > typemax(TN) || n < typemin(TN)) && continue
        data = open(readdlm, joinpath(@__DIR__, "data", "Sl$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            x        = row[1]
            expected = row[2]

            @test ClausenFunctions.sl(TN(n), x) == (-1)^TN(n)*ClausenFunctions.sl(TN(n), -x)
            @test ClausenFunctions.sl(TN(n), x) ≈ expected rtol=1e-14 atol=1e-14
            @test ClausenFunctions.sl(TN(n), Float16(x)) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
            @test ClausenFunctions.sl(TN(n), Float32(x)) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
            @test ClausenFunctions.sl(TN(n), Complex(x)) ≈ expected rtol=1e-14 atol=1e-14
            @test ClausenFunctions.sl(TN(n), Complex(Float16(x))) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
            @test ClausenFunctions.sl(TN(n), Complex(Float32(x))) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
        end
    end

    @test ClausenFunctions.sl(TN(1), 1//2) ≈ 1.3207963267948966 rtol=1e-14
    @test ClausenFunctions.sl(TN(2), 1//2) ≈ 0.92203590345077813 rtol=1e-14
    @test ClausenFunctions.sl(TN(3), 1//2) ≈ 0.63653415924141781 rtol=1e-14
    @test ClausenFunctions.sl(TN(4), 1//2) ≈ 0.90812931549667023 rtol=1e-14
    @test ClausenFunctions.sl(TN(5), 1//2) ≈ 0.51085256423059275 rtol=1e-14
    @test ClausenFunctions.sl(TN(6), 1//2) ≈ 0.88593812938731573 rtol=1e-14

    @test ClausenFunctions.sl(TN(1), Complex(1//2)) ≈ 1.3207963267948966 rtol=1e-14
    @test ClausenFunctions.sl(TN(2), Complex(1//2)) ≈ 0.92203590345077813 rtol=1e-14
    @test ClausenFunctions.sl(TN(3), Complex(1//2)) ≈ 0.63653415924141781 rtol=1e-14
    @test ClausenFunctions.sl(TN(4), Complex(1//2)) ≈ 0.90812931549667023 rtol=1e-14
    @test ClausenFunctions.sl(TN(5), Complex(1//2)) ≈ 0.51085256423059275 rtol=1e-14
    @test ClausenFunctions.sl(TN(6), Complex(1//2)) ≈ 0.88593812938731573 rtol=1e-14

    @test ClausenFunctions.sl(TN(1), 0.0 + 1.0im) ≈ 1.5707963267948966 - 0.5im rtol=1e-14
    @test ClausenFunctions.sl(TN(1), 0.0 - 1.0im) ≈ -1.5707963267948966 + 0.5im rtol=1e-14 # @todo check
    @test ClausenFunctions.sl(TN(2), 0.0 + 1.0im) ≈ 1.3949340668482264 - 1.5707963267948966im rtol=1e-14
    @test ClausenFunctions.sl(TN(2), 0.0 - 1.0im) ≈ 1.3949340668482264 - 1.5707963267948966im rtol=1e-14 # @todo check
    @test ClausenFunctions.sl(TN(2), 1.0 + 1.0im) ≈ 0.07413774005332982 - 1.07079632679489662im rtol=1e-14
    @test ClausenFunctions.sl(TN(2), 1.0 - 1.0im) ≈ 0.07413774005332982 + 1.07079632679489662im rtol=1e-14
    @test ClausenFunctions.sl(TN(3), 0.0 + 1.0im) ≈ 0.7853981633974483 + 1.5616007335148931im rtol=1e-14
    @test ClausenFunctions.sl(TN(3), 0.0 - 1.0im) ≈ -0.7853981633974483 - 1.5616007335148931im rtol=1e-14 # @todo check

    @test ClausenFunctions.sl(TN(-2), big(1))       == 0
    @test ClausenFunctions.sl(TN(-2), big(1) + 0im) == 0
    @test ClausenFunctions.sl(TN(-1), big(1))       == 0
    @test ClausenFunctions.sl(TN(-1), big(1) + 0im) == 0
    @test ClausenFunctions.sl(TN( 0), big(1))       ≈ -BigFloat("0.5") rtol=1e-40
    @test ClausenFunctions.sl(TN( 0), big(1) + 0im) ≈ -BigFloat("0.5") rtol=1e-40
    @test ClausenFunctions.sl(TN( 1), big(0))       == BigFloat(0)
    @test ClausenFunctions.sl(TN( 1), big(1))       ≈ BigFloat("1.07079632679489661923132169163975144209858") rtol=1e-40
    @test ClausenFunctions.sl(TN( 1), big(1) + 0im) ≈ BigFloat("1.07079632679489661923132169163975144209858") rtol=1e-40
    @test ClausenFunctions.sl(TN( 2), big(1))       ≈ BigFloat("0.324137740053329817241093475006273747120365") rtol=1e-40
    @test ClausenFunctions.sl(TN( 2), big(1) + 0im) ≈ BigFloat("0.324137740053329817241093475006273747120365") rtol=1e-40

    # test handling of negative zero
    for n in 1:2:101
        @test !signbit(ClausenFunctions.sl(TN(n), 0.0))
        @test signbit(ClausenFunctions.sl(TN(n), -0.0))
        @test !signbit(ClausenFunctions.sl(TN(n), BigFloat("0.0")))
        @test signbit(ClausenFunctions.sl(TN(n), BigFloat("-0.0")))

        @test !signbit(real(ClausenFunctions.sl(TN(n), Complex(0.0, 0.0))))
        @test !signbit(imag(ClausenFunctions.sl(TN(n), Complex(0.0, 0.0))))
        @test !signbit(real(ClausenFunctions.sl(TN(n), Complex(0.0, -0.0))))
        @test signbit(imag(ClausenFunctions.sl(TN(n), Complex(0.0, -0.0))))
        @test signbit(real(ClausenFunctions.sl(TN(n), Complex(-0.0, 0.0))))
        @test !signbit(imag(ClausenFunctions.sl(TN(n), Complex(-0.0, 0.0))))
        @test signbit(real(ClausenFunctions.sl(TN(n), Complex(-0.0, -0.0))))
        @test signbit(imag(ClausenFunctions.sl(TN(n), Complex(-0.0, -0.0))))
    end
end
