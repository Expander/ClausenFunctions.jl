@testset "TypeStability ($T)" for T in (Float16, Float32, Float64, BigFloat)
    @test @inferred(ClausenFunctions.cl1(T(1))) isa T
    @test @inferred(ClausenFunctions.cl2(T(1))) isa T
    
    if T != BigFloat
        @test @inferred(ClausenFunctions.cl3(T(1))) isa T
        @test @inferred(ClausenFunctions.cl4(T(1))) isa T
        @test @inferred(ClausenFunctions.cl5(T(1))) isa T
        @test @inferred(ClausenFunctions.cl6(T(1))) isa T
    end

    for n in -10:30
        @test @inferred(ClausenFunctions.cl(n,T(1))) isa T
        @test @inferred(ClausenFunctions.sl(n,T(1))) isa T
    end
end

@testset "TypeStability ($T)" for T in (ComplexF16, ComplexF32, ComplexF64, Complex{BigFloat})
    @test @inferred(ClausenFunctions.cl1(T(1))) isa T

    for n in -10:30
        @test @inferred(ClausenFunctions.cl(n,T(1))) isa T
        @test @inferred(ClausenFunctions.sl(n,T(1))) isa T
    end
end
