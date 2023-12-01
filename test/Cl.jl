@testset "cl" for TN in (Int8, Int16, Int32, Int64, Int128)
    for n in vcat(collect(1:16), [1000, 1001, 1_000_000, 0, -1, -2, -3])
        (n > typemax(TN) || n < typemin(TN)) && continue
        data = open(readdlm, joinpath(@__DIR__, "data", "Cl$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            x        = row[1]
            expected = row[2]
            sgn      = -(-1)^n

            cl  = ClausenFunctions.cl(TN(n), x)
            clm = sgn*ClausenFunctions.cl(TN(n), -x)
            ccl = ClausenFunctions.cl(TN(n), Complex(x))

            @test cl == clm
            @test cl ≈ expected rtol=1e-14 atol=1e-14
            @test ccl ≈ expected rtol=1e-14 atol=1e-14

            (n > typemax(TN) || n < typemin(TN)) && continue
            if n > 0
                @test ClausenFunctions.cl(TN(n), Float16(x)) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
                @test ClausenFunctions.cl(TN(n), Float32(x)) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
                @test ClausenFunctions.cl(TN(n), Complex(Float16(x))) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
                @test ClausenFunctions.cl(TN(n), Complex(Float32(x))) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
            end
        end

    end

    @test ClausenFunctions.cl(TN(1), 1//2) ≈ 0.70358563513784466 rtol=1e-14
    @test ClausenFunctions.cl(TN(2), 1//2) ≈ 0.84831187770367927 rtol=1e-14
    @test ClausenFunctions.cl(TN(3), 1//2) ≈ 0.92769631047023043 rtol=1e-14
    @test ClausenFunctions.cl(TN(4), 1//2) ≈ 0.54837172654589549 rtol=1e-14
    @test ClausenFunctions.cl(TN(5), 1//2) ≈ 0.89390286951083851 rtol=1e-14
    @test ClausenFunctions.cl(TN(6), 1//2) ≈ 0.49419627977618802 rtol=1e-14

    @test ClausenFunctions.cl(TN(1), Complex(1//2)) ≈ 0.70358563513784466 rtol=1e-14
    @test ClausenFunctions.cl(TN(2), Complex(1//2)) ≈ 0.84831187770367927 rtol=1e-14
    @test ClausenFunctions.cl(TN(3), Complex(1//2)) ≈ 0.92769631047023043 rtol=1e-14
    @test ClausenFunctions.cl(TN(4), Complex(1//2)) ≈ 0.54837172654589549 rtol=1e-14
    @test ClausenFunctions.cl(TN(5), Complex(1//2)) ≈ 0.89390286951083851 rtol=1e-14
    @test ClausenFunctions.cl(TN(6), Complex(1//2)) ≈ 0.49419627977618802 rtol=1e-14

    @test ClausenFunctions.cl(TN(-2), 1.0 + 1.0im) ≈ 0.50688457710655124 + 0.50950616274374456im rtol=1e-14
    @test ClausenFunctions.cl(TN(-2), 1.0 - 1.0im) ≈ 0.50688457710655124 - 0.50950616274374456im rtol=1e-14
    @test ClausenFunctions.cl(TN(-1), 1.0 + 1.0im) ≈ -0.08267495282794197 + 0.49171277760817954im rtol=1e-14
    @test ClausenFunctions.cl(TN(-1), 1.0 - 1.0im) ≈ -0.08267495282794197 - 0.49171277760817954im rtol=1e-14
    @test ClausenFunctions.cl(TN( 0), 1.0 + 1.0im) ≈ cot(1/2 + im/2)/2 rtol=1e-14
    @test ClausenFunctions.cl(TN( 0), 1.0 - 1.0im) ≈ coth(1/2 + im/2)*im/2 rtol=1e-14
    @test ClausenFunctions.cl(TN( 1), 0.0 + 1.0im) ≈ -0.0413248546129181 - 1.5707963267948966im rtol=1e-14
    @test ClausenFunctions.cl(TN( 1), 0.0 - 1.0im) ≈ -0.0413248546129181 - 1.5707963267948966im rtol=1e-14
    @test ClausenFunctions.cl(TN( 2), 0.0 + 1.0im) ≈ 1.5707963267948966 + 0.9861797794993302im rtol=1e-14
    @test ClausenFunctions.cl(TN( 2), 0.0 - 1.0im) ≈ -1.5707963267948966 - 0.9861797794993302im rtol=1e-14
    @test ClausenFunctions.cl(TN( 2), 1.0 + 1.0im) ≈ 1.4107754938116412 - 0.1044778629291566im rtol=1e-14
    @test ClausenFunctions.cl(TN( 2), 1.0 - 1.0im) ≈ 1.4107754938116412 + 0.1044778629291566im rtol=1e-14

    @test ClausenFunctions.cl(TN(-2), big(1))       ≈ 1im*(1 + exp(big(1)*im))*exp(big(1)*im)/(exp(big(1)*im) - 1)^3 rtol=1e-40
    @test ClausenFunctions.cl(TN(-2), big(1) + 0im) ≈ 1im*(1 + exp(big(1)*im))*exp(big(1)*im)/(exp(big(1)*im) - 1)^3 rtol=1e-40
    @test ClausenFunctions.cl(TN(-1), big(1))       ≈ exp(big(1)*im)/(exp(big(1)*im) - 1)^2 rtol=1e-40
    @test ClausenFunctions.cl(TN(-1), big(1) + 0im) ≈ exp(big(1)*im)/(exp(big(1)*im) - 1)^2 rtol=1e-40
    @test ClausenFunctions.cl(TN( 0), big(1))       ≈ cot(BigFloat("0.5"))/2 rtol=1e-40
    @test ClausenFunctions.cl(TN( 0), big(1) + 0im) ≈ cot(BigFloat("0.5"))/2 rtol=1e-40
    @test ClausenFunctions.cl(TN( 1), big(1))       ≈ BigFloat("0.04201950582536896172579838403790203712454") rtol=1e-40
    @test ClausenFunctions.cl(TN( 1), big(1) + 0im) ≈ BigFloat("0.04201950582536896172579838403790203712454") rtol=1e-40
    @test ClausenFunctions.cl(TN( 2), big(1))       ≈ BigFloat("1.0139591323607685042945743388859146875612") rtol=1e-40
    @test ClausenFunctions.cl(TN( 2), big(1) + 0im) ≈ BigFloat("1.0139591323607685042945743388859146875612") rtol=1e-40
end
