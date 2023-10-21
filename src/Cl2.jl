"""
    cl2(x::Real)::Real

Returns the value of the Clausen function ``\\operatorname{Cl}_2(x)``
for a real angle ``x`` of type `Real`.  The Clausen function is
defined as

```math
\\operatorname{Cl}_2(x) = -\\int_0^x \\log|2\\sin(t/2)| dt
```

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl2(1.0)
1.0139591323607684
```
"""
cl2(x::Real) = _cl2(float(x))

_cl2(x::Float16) = oftype(x, _cl2(Float32(x)))

_cl2(x::Float32) = oftype(x, _cl2(Float64(x)))

function _cl2(x::Float64)::Float64
    (x, sgn) = range_reduce_even(x)

    if iszero(x)
        zero(x)
    elseif x < pi/2
        P = (1.3888888888888889e-02, -4.3286930203743071e-04,
             3.2779814789973427e-06, -3.6001540369575084e-09)
        Q = (1.0000000000000000e+00, -3.6166589746694121e-02,
             3.6015827281202639e-04, -8.3646182842184428e-07)
        y = x*x
        y2 = y*y
        p = P[1] + y * P[2] + y2 * (P[3] + y * P[4])
        q = Q[1] + y * Q[2] + y2 * (Q[3] + y * Q[4])
        sgn*x*(1.0 - log(x) + y*p/q)
    elseif x < pi
        P = (6.4005702446195512e-01, -2.0641655351338783e-01,
             2.4175305223497718e-02, -1.2355955287855728e-03,
             2.5649833551291124e-05, -1.4783829128773320e-07)
        Q = (1.0000000000000000e+00, -2.5299102015666356e-01,
             2.2148751048467057e-02, -7.8183920462457496e-04,
             9.5432542196310670e-06, -1.8184302880448247e-08)
        y = pi - x
        z = y*y - pi*pi/8
        z2 = z*z
        z4 = z2*z2
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4]) + z4 * (P[5] + z * P[6])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4]) + z4 * (Q[5] + z * Q[6])
        sgn*y*p/q
    end
end
