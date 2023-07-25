@testset "cl1" begin
    data = open(readdlm, joinpath(@__DIR__, "data", "Cl1.txt"))

    for r in 1:size(data, 1)
        row      = data[r, :]
        x        = row[1]
        expected = row[2]

        @test ClausenFunctions.cl1(x) ≈ expected rtol=1e-13
        @test ClausenFunctions.cl1(-x) ≈ expected rtol=1e-13
        @test ClausenFunctions.cl1(x - 2pi) ≈ expected rtol=1e-12
        @test ClausenFunctions.cl1(x + 2pi) ≈ expected rtol=1e-13

        @test ClausenFunctions.cl1(Float16(x)) ≈ Float16(expected) atol=10*eps(Float16) rtol=30*eps(Float16)
        @test ClausenFunctions.cl1(Float32(x)) ≈ Float32(expected) atol=10*eps(Float32) rtol=30*eps(Float32)

        @test ClausenFunctions.cl1(x) ≈ real(ClausenFunctions.cl1(Complex(x))) rtol=1e-13
    end

    @test ClausenFunctions.cl1(1//2) ≈ 0.70358563513784466 rtol=1e-14
    @test ClausenFunctions.cl1(1) ≈ ClausenFunctions.cl1(1.0) rtol=1e-14
    @test ClausenFunctions.cl1(0.0) == Inf

    # test complex Clausen function cl1
    @test ClausenFunctions.cl1(Complex(1/2)) ≈ 0.70358563513784466 rtol=1e-14
    @test ClausenFunctions.cl1(Complex(1/2, 1)) ≈ -0.14296382186432422 - 1.06599896593689653im rtol=1e-14
    @test ClausenFunctions.cl1(Complex(1/2, -1)) ≈ -0.14296382186432422 + 1.06599896593689653im rtol=1e-14
    @test ClausenFunctions.cl1(Complex(1, 0)) ≈ 0.042019505825368962 rtol=1e-14
    @test ClausenFunctions.cl1(Complex(1, 1)) ≈ -0.34796082854253047 - 0.70210885509136190im rtol=1e-14
    @test ClausenFunctions.cl1(Complex(2, 0)) ≈ -0.52054343429085363 rtol=1e-14
    @test ClausenFunctions.cl1(Complex(1.0, 0.0)) ≈ 0.042019505825368962 rtol=1e-14
    @test ClausenFunctions.cl1(Complex(1.0, 1.0)) ≈ -0.34796082854253047 - 0.70210885509136190im rtol=1e-14
    @test ClausenFunctions.cl1(Complex(2.0, 0.0)) ≈ -0.52054343429085363 rtol=1e-14
    @test ClausenFunctions.cl1(Complex(2.0, 1.0)) ≈ -0.68284871442091651 - 0.28844676165196074im rtol=1e-14
    @test ClausenFunctions.cl1(Complex(2.0, -1.0)) ≈ -0.68284871442091651 + 0.28844676165196074im rtol=1e-14
    @test imag(ClausenFunctions.cl1(Complex(2.0, eps(Float64)))) < 0.0
    @test imag(ClausenFunctions.cl1(Complex(2.0, -eps(Float64)))) > 0.0
    @test ClausenFunctions.cl1(Complex(2.0, -eps(Float64))) ≈ -0.52054343429085363 rtol=1e-14
    @test ClausenFunctions.cl1(Complex(0)) == Inf
    @test ClausenFunctions.cl1(Complex(0.0)) == Inf
end
