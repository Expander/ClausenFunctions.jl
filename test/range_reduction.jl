@testset "range_reduction (n = $n)" for n in 1:2
    rr(x) = ClausenFunctions.range_reduce(n, x)[1]

    # cases with catastrophic cancellation (therefore not tested with rtol)
    @test rr(2pi)            ≈ zero(Float64)      atol=eps(Float64)
    @test rr(2pi + 1e-15)    ≈ 1e-15              atol=eps(Float64)
    @test rr(2pi + 1e-14)    ≈ 1e-14              atol=2*eps(Float64)
    @test rr(2pi + 1e-13)    ≈ 1e-13              atol=2*eps(Float64)
    @test rr(2pi + 1e-12)    ≈ 1e-12              atol=2*eps(Float64)

    # cases with catastrophic cancellation (therefore not tested with rtol)
    @test rr(prevfloat(2pi)) ≈ 1.4769252867665590e-15 atol=10*eps(Float64)
    @test rr(2pi - 1e-15)    ≈ 1e-15              atol=eps(Float64)
    @test rr(2pi - 1e-14)    ≈ 1e-14              atol=2*eps(Float64)
    @test rr(2pi - 1e-13)    ≈ 1e-13              atol=4*eps(Float64)
    @test rr(2pi - 1e-12)    ≈ 1e-12              atol=2*eps(Float64)

    # cases w/o catastrophic cancellation
    @test rr(3pi)            ≈ pi                 rtol=eps(Float64)
    @test rr(nextfloat(3pi)) ≈ 3.1415926535897920 rtol=eps(Float64)
    @test rr(3pi + 1e-15)    ≈ pi - 1e-15         rtol=eps(Float64)
    @test rr(3pi + 1e-14)    ≈ pi - 1e-14         rtol=eps(Float64)
    @test rr(3pi + 1e-13)    ≈ pi - 1e-13         rtol=2*eps(Float64)
    @test rr(3pi + 1e-12)    ≈ pi - 1e-12         rtol=eps(Float64)
end

@testset "two_pi_minus" begin
    f(x) = ClausenFunctions.two_pi_minus(x)

    # Float16
    @test f(Float16(2pi))               ≈ zero(Float16)        atol=2*eps(Float16)
    @test f(nextfloat(Float16(2pi)))    ≈ -eps(Float16(2pi))   atol=2*eps(Float16)
    @test f(nextfloat(Float16(2pi), 2)) ≈ -2*eps(Float16(2pi)) atol=2*eps(Float16)
    @test f(prevfloat(Float16(2pi)))    ≈ eps(Float16(2pi))    atol=2*eps(Float16)
    @test f(prevfloat(Float16(2pi), 2)) ≈ 2*eps(Float16(2pi))  atol=2*eps(Float16)

    # Float32
    @test f(Float32(2pi))               ≈ zero(Float32)        atol=2*eps(Float32)
    @test f(nextfloat(Float32(2pi)))    ≈ -eps(Float32(2pi))   atol=2*eps(Float32)
    @test f(nextfloat(Float32(2pi), 2)) ≈ -2*eps(Float32(2pi)) atol=2*eps(Float32)
    @test f(prevfloat(Float32(2pi)))    ≈ eps(Float32(2pi))    atol=2*eps(Float32)
    @test f(prevfloat(Float32(2pi), 2)) ≈ 2*eps(Float32(2pi))  atol=2*eps(Float32)

    # Float64
    @test f(2pi)               ≈ zero(Float64)    atol=2*eps(Float64)
    @test f(nextfloat(2pi))    ≈ -eps(2pi)        atol=2*eps(Float64)
    @test f(nextfloat(2pi, 2)) ≈ -2*eps(2pi)      atol=2*eps(Float64)
    @test f(prevfloat(2pi))    ≈ eps(2pi)         atol=2*eps(Float64)
    @test f(prevfloat(2pi, 2)) ≈ 2*eps(2pi)       atol=2*eps(Float64)
end
