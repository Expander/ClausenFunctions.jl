@testset "range_reduction_odd" begin
    rr(n, x) = ClausenFunctions.range_reduce(n, x)[1]

    # cases with catastrophic cancellation (therefore not tested with rtol)
    @test rr(1, 2pi)            ≈ zero(Float64)      atol=eps(Float64)
    @test rr(1, 2pi + 1e-15)    ≈ 1e-15              atol=eps(Float64)
    @test rr(1, 2pi + 1e-14)    ≈ 1e-14              atol=2*eps(Float64)
    @test rr(1, 2pi + 1e-13)    ≈ 1e-13              atol=2*eps(Float64)
    @test rr(1, 2pi + 1e-12)    ≈ 1e-12              atol=2*eps(Float64)

    # cases with catastrophic cancellation (therefore not tested with rtol)
    @test rr(1, prevfloat(2pi)) ≈ 1.4769252867665590e-15 atol=10*eps(Float64)
    @test rr(1, 2pi - 1e-15)    ≈ 1e-15              atol=eps(Float64)
    @test rr(1, 2pi - 1e-14)    ≈ 1e-14              atol=2*eps(Float64)
    @test rr(1, 2pi - 1e-13)    ≈ 1e-13              atol=4*eps(Float64)
    @test rr(1, 2pi - 1e-12)    ≈ 1e-12              atol=2*eps(Float64)

    # cases w/o catastrophic cancellation
    @test rr(1, 3pi)            ≈ pi                 rtol=eps(Float64)
    @test rr(1, nextfloat(3pi)) ≈ 3.1415926535897920 rtol=eps(Float64)
    @test rr(1, 3pi + 1e-15)    ≈ pi - 1e-15         rtol=eps(Float64)
    @test rr(1, 3pi + 1e-14)    ≈ pi - 1e-14         rtol=eps(Float64)
    @test rr(1, 3pi + 1e-13)    ≈ pi - 1e-13         rtol=2*eps(Float64)
    @test rr(1, 3pi + 1e-12)    ≈ pi - 1e-12         rtol=eps(Float64)
end
