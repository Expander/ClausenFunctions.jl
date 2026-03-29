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
    ref_256 = setprecision(BigFloat, 256) do
        ClausenFunctions.cl2(BigFloat(pi)/4)
    end
    ref_512 = setprecision(BigFloat, 512) do
        ClausenFunctions.cl2(BigFloat(pi)/4)
    end

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

    for (r256, r512) in results
        @test BigFloat(r256, precision=256) == BigFloat(ref_256, precision=256)
        @test BigFloat(r512, precision=512) == BigFloat(ref_512, precision=512)
    end
end
