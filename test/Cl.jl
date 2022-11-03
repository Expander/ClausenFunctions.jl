@testset "cl" begin
    for n in vcat(collect(1:16), [1000, 1001, 1_000_000])
        data = open(readdlm, joinpath(@__DIR__, "data", "Cl$(n).txt"))

        for r in 1:size(data, 1)
            row      = data[r, :]
            x        = row[1]
            expected = row[2]
            sgn      = -(-1)^n

            cl  = ClausenFunctions.cl(n, x)
            clm = sgn*ClausenFunctions.cl(n, -x)

            @test cl == clm
            @test cl ≈ expected atol=1e-14
            @test ClausenFunctions.cl(n, Float16(x)) ≈ Float16(expected) atol=30*eps(Float16) rtol=30*eps(Float16)
            @test ClausenFunctions.cl(n, Float32(x)) ≈ Float32(expected) atol=30*eps(Float32) rtol=30*eps(Float32)
        end

    end

    @test ClausenFunctions.cl(1, 1//2) ≈ 0.70358563513784466 atol=1e-14
    @test ClausenFunctions.cl(2, 1//2) ≈ 0.84831187770367927 atol=1e-14
    @test ClausenFunctions.cl(3, 1//2) ≈ 0.92769631047023043 atol=1e-14
    @test ClausenFunctions.cl(4, 1//2) ≈ 0.54837172654589549 atol=1e-14
    @test ClausenFunctions.cl(5, 1//2) ≈ 0.89390286951083851 atol=1e-14
    @test ClausenFunctions.cl(6, 1//2) ≈ 0.49419627977618802 atol=1e-14

    @test_throws DomainError ClausenFunctions.cl(0, 1.0)
    @test_throws DomainError ClausenFunctions.cl(-1, 1.0)
    @test_throws DomainError ClausenFunctions.cl(-2, 1.0)
end
