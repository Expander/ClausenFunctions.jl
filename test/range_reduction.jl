@testset "range_reduction (n = $n)" for n in 1:2
    rr(x) = ClausenFunctions.range_reduce(n, x)[1]

    # cases with catastrophic cancellation (therefore not tested with rtol)
    @test rr(2pi)            ≈ zero(Float64)      atol=eps(Float64)
    @test rr(2pi + 1e-15)    ≈ 1e-15              atol=eps(Float64)
    @test rr(2pi + 1e-14)    ≈ 1e-14              atol=2*eps(Float64)
    @test rr(2pi + 1e-13)    ≈ 1e-13              atol=2*eps(Float64)
    @test rr(2pi + 1e-12)    ≈ 1e-12              atol=2*eps(Float64)

    # cases with catastrophic cancellation (therefore not tested with rtol)
    @test rr(prevfloat(2pi)) ≈ 1.4769252867665590e-15 atol=10*eps(Float64)
    @test rr(2pi - 1e-15)    ≈ 1e-15              atol=eps(Float64)
    @test rr(2pi - 1e-14)    ≈ 1e-14              atol=2*eps(Float64)
    @test rr(2pi - 1e-13)    ≈ 1e-13              atol=4*eps(Float64)
    @test rr(2pi - 1e-12)    ≈ 1e-12              atol=2*eps(Float64)

    # cases w/o catastrophic cancellation
    @test rr(3pi)            ≈ pi                 rtol=eps(Float64)
    @test rr(nextfloat(3pi)) ≈ 3.1415926535897920 rtol=eps(Float64)
    @test rr(3pi + 1e-15)    ≈ pi - 1e-15         rtol=eps(Float64)
    @test rr(3pi + 1e-14)    ≈ pi - 1e-14         rtol=eps(Float64)
    @test rr(3pi + 1e-13)    ≈ pi - 1e-13         rtol=2*eps(Float64)
    @test rr(3pi + 1e-12)    ≈ pi - 1e-12         rtol=eps(Float64)
end

@testset "range_reduction BigFloat (n = $n)" for n in 1:2
    rr(x) = ClausenFunctions.range_reduce(n, x)

    setprecision(BigFloat, 256) do
        # reference: reduce x to [-π, π] at much higher precision
        function ref(x)
            r = setprecision(BigFloat, 4096) do
                y = mod(x, 2*BigFloat(pi))
                y > BigFloat(pi) ? y - 2*BigFloat(pi) : y
            end
            BigFloat(abs(r)), (iseven(n) && signbit(r)) ? -1 : 1
        end

        twopi = 2*BigFloat(pi)
        for x in (big"0.5", big"-0.5", big"3.0", big"-3.0", big"4.0", big"-4.0",
                  big"1e100", -big"1e30", big"123456.789",
                  twopi*1000 + big"1e-40", twopi - big"1e-50",
                  -(twopi - big"1e-50"), 7*BigFloat(pi) + big"1e-60")
            (y, sgn) = rr(x)
            (yr, sgnr) = ref(x)
            @test y ≈ yr rtol=4*eps(BigFloat)
            @test sgn == sgnr
            @test 0 <= y <= BigFloat(pi)
        end

        # zero, infinity and NaN
        @test rr(big"0.0") == (0, 1)
        @test isnan(rr(BigFloat(Inf))[1])
        @test isnan(rr(BigFloat(-Inf))[1])
        @test isnan(rr(BigFloat(NaN))[1])
    end
end

@testset "two_pi_minus" begin
    f(x) = ClausenFunctions.two_pi_minus(x)

    # Float16
    @test f(Float16(2pi))               ≈ zero(Float16)        atol=2*eps(Float16)
    @test f(nextfloat(Float16(2pi)))    ≈ -eps(Float16(2pi))   atol=2*eps(Float16)
    @test f(nextfloat(Float16(2pi), 2)) ≈ -2*eps(Float16(2pi)) atol=2*eps(Float16)
    @test f(prevfloat(Float16(2pi)))    ≈ eps(Float16(2pi))    atol=2*eps(Float16)
    @test f(prevfloat(Float16(2pi), 2)) ≈ 2*eps(Float16(2pi))  atol=2*eps(Float16)

    # Float32
    @test f(Float32(2pi))               ≈ zero(Float32)        atol=2*eps(Float32)
    @test f(nextfloat(Float32(2pi)))    ≈ -eps(Float32(2pi))   atol=2*eps(Float32)
    @test f(nextfloat(Float32(2pi), 2)) ≈ -2*eps(Float32(2pi)) atol=2*eps(Float32)
    @test f(prevfloat(Float32(2pi)))    ≈ eps(Float32(2pi))    atol=2*eps(Float32)
    @test f(prevfloat(Float32(2pi), 2)) ≈ 2*eps(Float32(2pi))  atol=2*eps(Float32)

    # Float64
    @test f(2pi)               ≈ zero(Float64)    atol=2*eps(Float64)
    @test f(nextfloat(2pi))    ≈ -eps(2pi)        atol=2*eps(Float64)
    @test f(nextfloat(2pi, 2)) ≈ -2*eps(2pi)      atol=2*eps(Float64)
    @test f(prevfloat(2pi))    ≈ eps(2pi)         atol=2*eps(Float64)
    @test f(prevfloat(2pi, 2)) ≈ 2*eps(2pi)       atol=2*eps(Float64)

    # BigFloat
    twopi = 2*BigFloat(pi)
    @test f(twopi)                       ≈ zero(BigFloat) atol=2*eps(BigFloat)
    @test f(nextfloat(twopi))            ≈ -eps(twopi)    atol=2*eps(BigFloat)
    @test f(nextfloat(nextfloat(twopi))) ≈ -2*eps(twopi)  atol=2*eps(BigFloat)
    @test f(prevfloat(twopi))            ≈ eps(twopi)     atol=2*eps(BigFloat)
    @test f(prevfloat(prevfloat(twopi))) ≈ 2*eps(twopi)   atol=2*eps(BigFloat)
end

@testset "rem_pio2_large" begin
    # reference: x - k π/2 with k = round(x 2/π), and k mod 4
    function ref(x)
        setprecision(BigFloat, 3000) do
            X = BigFloat(x)
            k = round(X/(BigFloat(pi)/2))
            (Int(mod(k, 4)), X - k*BigFloat(pi)/2)
        end
    end
    # error of the angle n π/2 + y, modulo 2π, relative to the reduced argument
    function relerr(x)
        (n, yh, yl) = ClausenFunctions.rem_pio2_large(x)
        (nr, yr) = ref(x)
        setprecision(BigFloat, 3000) do
            d = (n - nr)*BigFloat(pi)/2 + (BigFloat(yh) + BigFloat(yl)) - yr
            d -= round(d/(2*BigFloat(pi)))*2*BigFloat(pi)
            Float64(abs(d)/abs(yr))
        end
    end

    xs = [s*ldexp(1 + k/8, e) for e in 14:1023 for k in 0:7 for s in (-1, 1)]
    append!(xs, [2.0^14, prevfloat(2.0^15), floatmax(Float64), -floatmax(Float64),
                 6381956970095103*2.0^797, -6381956970095103*2.0^797]) # closest to a multiple of π/2
    for x in xs
        @test relerr(x) < 2.0^-60
        (n, yh, yl) = ClausenFunctions.rem_pio2_large(x)
        @test abs(yh) <= 0.786
        @test abs(yl) <= eps(yh)/2
    end
end

@testset "Payne-Hanek table" begin
    # _PH_TABLE is built at load (precompile) time; rebuild it at run time
    @test ClausenFunctions._ph_table() == ClausenFunctions._PH_TABLE
    # D₀, D₁, D₂ have at most 51 significant bits (the last two of the
    # 52 stored significand bits are zero)
    for (D0, D1, D2, D3) in ClausenFunctions._PH_TABLE, D in (D0, D1, D2)
        @test reinterpret(UInt64, D) & 0x3 == 0
    end
end
