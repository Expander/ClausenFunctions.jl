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
    # Reference values computed with PolyLog.li(2, cis(θ)) at 1024-bit precision.

    setprecision(BigFloat, 512) do
        @test ClausenFunctions.cl2(BigFloat(pi)/6)    ≈ big"0.86437913105389274962503631519021947866818857640368970418203768977532471558206" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,BigFloat(pi)/6)   ≈ big"0.86437913105389274962503631519021947866818857640368970418203768977532471558206" rtol=big"1e-60"
        @test ClausenFunctions.cl2(BigFloat(pi)/4)    ≈ big"0.98187215105020335671792454306019566713079097166073046157661313465315566504977" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,BigFloat(pi)/4)   ≈ big"0.98187215105020335671792454306019566713079097166073046157661313465315566504977" rtol=big"1e-60"
        @test ClausenFunctions.cl2(BigFloat(pi)/3)    ≈ big"1.01494160640965362502120255427452028594168930753029979201748910677659747625824" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,BigFloat(pi)/3)   ≈ big"1.01494160640965362502120255427452028594168930753029979201748910677659747625824" rtol=big"1e-60"
        @test ClausenFunctions.cl2(BigFloat(pi)/2)    ≈ big"0.91596559417721901505460351493238411077414937428167213426649811962176301977625" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,BigFloat(pi)/2)   ≈ big"0.91596559417721901505460351493238411077414937428167213426649811962176301977625" rtol=big"1e-60"
        @test ClausenFunctions.cl2(big"1")            ≈ big"1.01395913236076850429457433888591468756117928007771731687704851226813781234608" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,big"1")           ≈ big"1.01395913236076850429457433888591468756117928007771731687704851226813781234608" rtol=big"1e-60"
        @test ClausenFunctions.cl2(2*BigFloat(pi)/3)  ≈ big"0.67662773760643575001413503618301352396112620502019986134499273785106498417216" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,2*BigFloat(pi)/3) ≈ big"0.67662773760643575001413503618301352396112620502019986134499273785106498417216" rtol=big"1e-60"
        @test ClausenFunctions.cl2(big"0.0")  == big"0.0"
        @test ClausenFunctions.cl(2,big"0.0") == big"0.0"
        @test ClausenFunctions.cl2(BigFloat(pi))  ≈ big"0.0" atol=big"1e-60"
        @test ClausenFunctions.cl(2,BigFloat(pi)) ≈ big"0.0" atol=big"1e-60"

        # Antisymmetry: Cl₂(-x) = -Cl₂(x)
        @test ClausenFunctions.cl2(-BigFloat(pi)/4)  ≈ -ClausenFunctions.cl2(BigFloat(pi)/4) rtol=big"1e-60"
        @test ClausenFunctions.cl(2,-BigFloat(pi)/4) ≈ -ClausenFunctions.cl2(BigFloat(pi)/4) rtol=big"1e-60"

        # Periodicity: Cl₂(x + 2pi) = Cl₂(x)
        @test ClausenFunctions.cl2(big"1" + 2*BigFloat(pi))  ≈ ClausenFunctions.cl2(big"1") rtol=big"1e-40"
        @test ClausenFunctions.cl(2,big"1" + 2*BigFloat(pi)) ≈ ClausenFunctions.cl2(big"1") rtol=big"1e-40"
    end

    # Test at higher precision (1000 bits ≈ 300 digits)
    setprecision(BigFloat, 1000) do
        @test ClausenFunctions.cl2(BigFloat(pi)/2)  ≈ big"0.91596559417721901505460351493238411077414937428167213426649811962176301977625" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,BigFloat(pi)/2) ≈ big"0.91596559417721901505460351493238411077414937428167213426649811962176301977625" rtol=big"1e-60"
        @test ClausenFunctions.cl2(big"1")          ≈ big"1.01395913236076850429457433888591468756117928007771731687704851226813781234608" rtol=big"1e-60"
        @test ClausenFunctions.cl(2,big"1")         ≈ big"1.01395913236076850429457433888591468756117928007771731687704851226813781234608" rtol=big"1e-60"
    end
end
