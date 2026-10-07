if isdefined(Base, :get_extension)

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

end
