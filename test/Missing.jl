@testset "Missing" begin
    @test ismissing(ClausenFunctions.cl1(missing))
    @test ismissing(ClausenFunctions.cl2(missing))
    @test ismissing(ClausenFunctions.cl3(missing))
    @test ismissing(ClausenFunctions.cl4(missing))
    @test ismissing(ClausenFunctions.cl5(missing))
    @test ismissing(ClausenFunctions.cl6(missing))

    for n in -16:16
        @test ismissing(ClausenFunctions.cl(missing, 1.0))
        @test ismissing(ClausenFunctions.cl(n, missing))
        @test ismissing(ClausenFunctions.sl(missing, 1.0))
        @test ismissing(ClausenFunctions.sl(n, missing))
    end
end
