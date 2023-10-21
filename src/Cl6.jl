"""
    cl6(x::Real)::Real

Returns the value of the Clausen function ``\\operatorname{Cl}_6(x)``
for a real angle ``x`` of type `Real`.  This function is defined as

```math
\\operatorname{Cl}_6(x) = \\Im[\\operatorname{Li}_6(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\sin(kx)}{k^6}
```

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl6(1.0)
0.855629273183937
```
"""
cl6(x::Real) = _cl6(float(x))

_cl6(x::Float16) = oftype(x, _cl6(Float32(x)))

_cl6(x::Float32) = oftype(x, _cl6(Float64(x)))

function _cl6(x::Float64)::Float64
    (x, sgn) = range_reduce_even(x)

    if iszero(x)
        zero(x)
    elseif x < pi/2
        P = (1.0369277551433699e+00, -2.087195444107175e-01,
             2.0652251045312954e-02, -1.383438138256840e-04)
        Q = (1.0000000000000000e+00, -8.0784096827362542e-03,
             5.8074568862993102e-06, -5.1960620033050114e-10)
        y = x*x
        y2 = y*y
        p = P[1] + y * P[2] + y2 * (P[3] + y * P[4])
        q = Q[1] + y * Q[2] + y2 * (Q[3] + y * Q[4])
        sgn*x*(p/q - 1/120*y2*log(x))
    else
        P = (7.9544504578027050e-01, -1.9255025309738589e-01,
             1.5805208288846591e-02, -5.4175380521534706e-04,
             6.7577493541009068e-06)
        Q = (1.0000000000000000e+00, -7.0798422394109274e-02,
             7.1744189715634762e-04,  3.9098747334347093e-06,
             3.5669441618295266e-08,  2.5315391843409925e-10)
        y = pi - x
        z = y*y - pi*pi/8
        z2 = z*z
        z4 = z2*z2
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4]) + z4 * P[5]
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4]) + z4 * (Q[5] + z * Q[6])
        sgn*y*p/q
    end
end
