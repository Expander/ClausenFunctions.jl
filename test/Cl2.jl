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
    # Float64(π) = π - 1.2246467991473532e-16, so Cl₂(Float64(π)) ≈ 1.2246e-16 log(2)
    @test ClausenFunctions.cl2(pi) ≈ 8.488604760107494e-17 rtol=2*eps(Float64)
    @test ClausenFunctions.cl2(-pi) ≈ -8.488604760107494e-17 rtol=2*eps(Float64)
    @test ClausenFunctions.cl2(pi/2) ≈ 0.915965594177219015054603514932384110 rtol=1e-14
    @test ClausenFunctions.cl2(1//2) ≈ 0.84831187770367927 rtol=1e-14

    # test handling of negative zero
    @test !signbit(ClausenFunctions.cl2(0.0))
    @test signbit(ClausenFunctions.cl2(-0.0))
end



@testset "cl2 Float64 accuracy" begin
    # reference value of Cl₂ for the exact Float64 argument x
    ref(x) = setprecision(BigFloat, 256) do
        Float64(ClausenFunctions.cl2(BigFloat(x)))
    end
    # error in units of the last place of the reference value
    ulps(x) = (r = ref(x); abs(ClausenFunctions.cl2(x) - r)/eps(r))

    # random arguments, also beyond [-π,π]
    for x in range(-4pi, stop=4pi, length=2001)
        @test ulps(x) < 2
    end

    # arguments near 0, π, 2π, 3π and large multiples of π, approached
    # from both sides: the result has a small magnitude there, so
    # argument reduction errors would be amplified
    for k in (0, 1, 2, 3, -1, -2, 101, 2^20, 10^8), e in -15:-1, s in (-1, 1)
        x = k*pi + s*10.0^e
        iszero(x) && continue
        @test ulps(x) < 2
    end

    # large arguments (reference with enough bits to reduce |x| ≤ 1e300)
    ref_big(x) = setprecision(BigFloat, 2048) do
        Float64(ClausenFunctions.cl2(BigFloat(x)))
    end
    for x in (1e3, -1e3, 12345.678, 1e6, 1e8, -1e15, 1e100, 1e300)
        r = ref_big(x)
        @test abs(ClausenFunctions.cl2(x) - r)/eps(r) < 2
    end

    # tiny and subnormal arguments
    for x in (1e-300, 5e-324, nextfloat(0.0, 100), 2.0^-1000)
        @test ulps(x) < 2
        @test ClausenFunctions.cl2(-x) == -ClausenFunctions.cl2(x)
    end

    @test isnan(ClausenFunctions.cl2(NaN))
    @test isnan(ClausenFunctions.cl2(Inf))
    @test isnan(ClausenFunctions.cl2(-Inf))
end


@testset "cl2 BigFloat large arguments" begin
    # Cl₂ is 2π-periodic; the argument reduction must not lose precision
    # for large |x| or near multiples of 2π
    setprecision(BigFloat, 256) do
        twopi = 2*BigFloat(pi)
        for x in (big"1e100", -big"1e30", big"123456.789", twopi*1000 + big"1e-40",
                  twopi - big"1e-50", 7*BigFloat(pi) + big"1e-60")
            # reference: reduce the exact x to [0, 2π) at much higher precision
            r = setprecision(BigFloat, 4096) do
                mod(BigFloat(x; precision=4096), 2*BigFloat(pi))
            end
            @test ClausenFunctions.cl2(x) ≈ ClausenFunctions.cl2(r) rtol=10*eps(BigFloat)
        end
        # input with more precision than the working precision
        x = setprecision(BigFloat, 2000) do
            2*BigFloat(pi)*1000 + big"1e-200"
        end
        @test ClausenFunctions.cl2(x) ≈ ClausenFunctions.cl2(big"1e-200") rtol=10*eps(BigFloat)
    end
end


@testset "cl2 ForwardDiff" begin
    # d/dx Cl₂(x) = Cl₁(x) = -log|2 sin(x/2)|
    for x in (0.1, 1.0, 2.0, 3.0, -1.0, 7.0, 1e3)
        @test ForwardDiff.derivative(ClausenFunctions.cl2, x) ≈ ClausenFunctions.cl1(x) rtol=1e-14
        @test ForwardDiff.derivative(ClausenFunctions.cl2, x) ≈ -log(abs(2*sin(x/2))) rtol=1e-12
    end
    # second derivative: d²/dx² Cl₂(x) = -cot(x/2)/2
    for x in (0.1, 1.0, 2.0, 3.0)
        d2 = ForwardDiff.derivative(y -> ForwardDiff.derivative(ClausenFunctions.cl2, y), x)
        @test d2 ≈ -cot(x/2)/2 rtol=1e-12
    end
    # BigFloat dual numbers
    setprecision(BigFloat, 256) do
        x = BigFloat(1)
        @test ForwardDiff.derivative(ClausenFunctions.cl2, x) ≈ ClausenFunctions.cl1(x) rtol=10*eps(BigFloat)
    end
    # Float32
    @test ForwardDiff.derivative(ClausenFunctions.cl2, 1.0f0) isa Float32
end


@testset "cl2 BigFloat" begin
    decimal_digits = 100
    binary_digits = ceil(Int, decimal_digits*log(10)/log(2))

    setprecision(BigFloat, binary_digits) do
        eps = 10.0^(1 - decimal_digits)

        @test ClausenFunctions.cl2(BigFloat(pi)/6)    ≈ BigFloat("0.8643791310538927496250363151902194786681885764036897041820376897753247155820641851870219305078075779") rtol=eps
        @test ClausenFunctions.cl(2,BigFloat(pi)/6)   ≈ BigFloat("0.8643791310538927496250363151902194786681885764036897041820376897753247155820641851870219305078075779") rtol=eps
        @test ClausenFunctions.cl2(BigFloat(pi)/4)    ≈ BigFloat("0.9818721510502033567179245430601956671307909716607304615766131346531556650497696362249028028843877241") rtol=eps
        @test ClausenFunctions.cl(2,BigFloat(pi)/4)   ≈ BigFloat("0.9818721510502033567179245430601956671307909716607304615766131346531556650497696362249028028843877241") rtol=eps
        @test ClausenFunctions.cl2(BigFloat(pi)/3)    ≈ BigFloat("1.0149416064096536250212025542745202859416893075302997920174891067765974762582440221364703542282566949") rtol=eps
        @test ClausenFunctions.cl(2,BigFloat(pi)/3)   ≈ BigFloat("1.0149416064096536250212025542745202859416893075302997920174891067765974762582440221364703542282566949") rtol=eps
        @test ClausenFunctions.cl2(BigFloat(pi)/2)    ≈ BigFloat("0.9159655941772190150546035149323841107741493742816721342664981196217630197762547694793565129261151062") rtol=eps
        @test ClausenFunctions.cl(2,BigFloat(pi)/2)   ≈ BigFloat("0.9159655941772190150546035149323841107741493742816721342664981196217630197762547694793565129261151062") rtol=eps
        @test ClausenFunctions.cl2(BigFloat("1"))     ≈ BigFloat("1.0139591323607685042945743388859146875611792800777173168770485122681378123460795573363882186547712204") rtol=eps
        @test ClausenFunctions.cl(2,BigFloat("1"))    ≈ BigFloat("1.0139591323607685042945743388859146875611792800777173168770485122681378123460795573363882186547712204") rtol=eps
        @test ClausenFunctions.cl2(2*BigFloat(pi)/3)  ≈ BigFloat("0.6766277376064357500141350361830135239611262050201998613449927378510649841721626814243135694855044633") rtol=eps
        @test ClausenFunctions.cl(2,2*BigFloat(pi)/3) ≈ BigFloat("0.6766277376064357500141350361830135239611262050201998613449927378510649841721626814243135694855044633") rtol=eps
        @test ClausenFunctions.cl2(BigFloat("0.0"))  == BigFloat("0.0")
        @test ClausenFunctions.cl(2,BigFloat("0.0")) == BigFloat("0.0")
        @test ClausenFunctions.cl2(BigFloat(pi))  ≈ BigFloat("0.0") atol=eps
        @test ClausenFunctions.cl(2,BigFloat(pi)) ≈ BigFloat("0.0") atol=eps

        # Antisymmetry: Cl₂(-x) = -Cl₂(x)
        @test ClausenFunctions.cl2(-BigFloat(pi)/4)  ≈ -ClausenFunctions.cl2(BigFloat(pi)/4) rtol=eps
        @test ClausenFunctions.cl(2,-BigFloat(pi)/4) ≈ -ClausenFunctions.cl2(BigFloat(pi)/4) rtol=eps

        # Periodicity: Cl₂(x + 2pi) = Cl₂(x)
        @test ClausenFunctions.cl2(BigFloat("1") + 2*BigFloat(pi))  ≈ ClausenFunctions.cl2(BigFloat("1")) rtol=eps
        @test ClausenFunctions.cl(2,BigFloat("1") + 2*BigFloat(pi)) ≈ ClausenFunctions.cl2(BigFloat("1")) rtol=eps
    end

    # Test at higher precision
    decimal_digits = 200
    binary_digits = ceil(Int, decimal_digits*log(10)/log(2))

    setprecision(BigFloat, binary_digits) do
        eps = 10.0^(1 - decimal_digits)

        @test ClausenFunctions.cl2(BigFloat(pi)/2)  ≈ BigFloat("0.91596559417721901505460351493238411077414937428167213426649811962176301977625476947935651292611510624857442261919619957903589880332585905943159473748115840699533202877331946051903872747816408786590902") rtol=eps
        @test ClausenFunctions.cl(2,BigFloat(pi)/2) ≈ BigFloat("0.91596559417721901505460351493238411077414937428167213426649811962176301977625476947935651292611510624857442261919619957903589880332585905943159473748115840699533202877331946051903872747816408786590902") rtol=eps
        @test ClausenFunctions.cl2(BigFloat("1"))   ≈ BigFloat("1.01395913236076850429457433888591468756117928007771731687704851226813781234607955733638821865477122042157440086434150311308425232269856980238591034316884447292797795707328765622668664434914352893878283") rtol=eps
        @test ClausenFunctions.cl(2,BigFloat("1"))  ≈ BigFloat("1.01395913236076850429457433888591468756117928007771731687704851226813781234607955733638821865477122042157440086434150311308425232269856980238591034316884447292797795707328765622668664434914352893878283") rtol=eps
    end
end


@testset "cl2 BigFloat thread safety" begin
    # Hammer _cl2_ensure_coeffs from multiple threads with varying
    # precisions to exercise the cache lock.

    # reference value, calculated sequentially in the main thread
    ref_256 = setprecision(BigFloat, 256) do
        ClausenFunctions.cl2(BigFloat(pi)/4)
    end
    ref_512 = setprecision(BigFloat, 512) do
        ClausenFunctions.cl2(BigFloat(pi)/4)
    end

    # values calculated in parallel in multiple threads
    results = Vector{Tuple{BigFloat,BigFloat}}(undef, 64)
    Threads.@threads for i in 1:64
        r256 = setprecision(BigFloat, 256) do
            ClausenFunctions.cl2(BigFloat(pi)/4)
        end
        r512 = setprecision(BigFloat, 512) do
            ClausenFunctions.cl2(BigFloat(pi)/4)
        end
        results[i] = (r256, r512)
    end

    # test equality
    for (r256, r512) in results
        setprecision(BigFloat, 256) do
            @test BigFloat(r256) == BigFloat(ref_256)
        end
        setprecision(BigFloat, 512) do
            @test BigFloat(r512) == BigFloat(ref_512)
        end
    end
end
