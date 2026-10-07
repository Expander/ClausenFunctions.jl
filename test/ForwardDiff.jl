if isdefined(Base, :get_extension)

    cl0(x) = cot(x/2)/2

    @testset "cl1 ForwardDiff ($T)" for T in (Float16, Float32, Float64, BigFloat)
        for x in (0.1, 1.0, 2.0, 3.0, -1.0, 7.0, 1e3)
            @test ForwardDiff.derivative(ClausenFunctions.cl1, T(x)) ≈ -cl0(T(x)) rtol=eps(T)
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
            @test d2 ≈ -cl0(x) rtol=1e-12
        end
        # BigFloat dual numbers
        setprecision(BigFloat, 256) do
            x = BigFloat(1)
            @test ForwardDiff.derivative(ClausenFunctions.cl2, x) ≈ ClausenFunctions.cl1(x) rtol=10*eps(BigFloat)
        end
        # Float32
        @test ForwardDiff.derivative(ClausenFunctions.cl2, 1.0f0) isa Float32
    end

    @testset "cl3 ForwardDiff (Float64)" begin
        for x in (0.1, 1.0, 2.0, 3.0, -1.0, 7.0, 1e3)
            @test ForwardDiff.derivative(ClausenFunctions.cl3, x) ≈ -ClausenFunctions.cl2(x) rtol=eps(Float64)
        end
    end

    @testset "cl4 ForwardDiff (Float64)" begin
        for x in (0.1, 1.0, 2.0, 3.0, -1.0, 7.0, 1e3)
            @test ForwardDiff.derivative(ClausenFunctions.cl4, x) ≈ ClausenFunctions.cl3(x) rtol=eps(Float64)
        end
    end

    @testset "cl5 ForwardDiff (Float64)" begin
        for x in (0.1, 1.0, 2.0, 3.0, -1.0, 7.0, 1e3)
            @test ForwardDiff.derivative(ClausenFunctions.cl5, x) ≈ -ClausenFunctions.cl4(x) rtol=eps(Float64)
        end
    end

    @testset "cl6 ForwardDiff (Float64)" begin
        for x in (0.1, 1.0, 2.0, 3.0, -1.0, 7.0, 1e3)
            @test ForwardDiff.derivative(ClausenFunctions.cl6, x) ≈ ClausenFunctions.cl5(x) rtol=eps(Float64)
        end
    end

end
