"""
    cl4(x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Cl}_4(x)``
for a real angle ``x`` of type `Float64`.  This function is defined as

```math
\\operatorname{Cl}_4(x) = \\Im[\\operatorname{Li}_4(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\sin(kx)}{k^4}
```

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl4(1.0)
0.8958052386793799
```
"""
function cl4(x::Float64)::Float64
    zeta3 = 1.2020569031595943
    pi28 = pi*pi/8.0

    (x, sgn) = range_reduce_even(x)

    if x == zero(x) || x == pi
        return zero(x)
    end

    if x < 0.5*pi
        P = (-3.0555555555555556e-01,  6.0521392328447206e-03,
             -1.9587493942041528e-05, -3.1137343767030358e-08)
        Q = (1.0000000000000000e+00, -2.2079728398400851e-02,
             1.0887447112236682e-04, -6.1847621370547954e-08)
        y = x*x
        y2 = y*y
        p = P[1] + y * P[2] + y2 * (P[3] + y * P[4])
        q = Q[1] + y * Q[2] + y2 * (Q[3] + y * Q[4])
        sgn*x*(zeta3 + y*(p/q + 1/6*log(x)))
    else
        P = (7.6223911686491336e-01, -2.4339587368267260e-01,
             2.8715364937979943e-02, -1.5368612510964667e-03,
             3.6261044225761673e-05, -2.8557977333851308e-07)
        Q = (1.0000000000000000e+00, -1.7465715261403233e-01,
             9.5439417991615653e-03, -1.7325070821666274e-04,
             5.9283675098376635e-07,  9.4127575773361230e-10)
        y = pi - x
        z = y*y - pi28
        z2 = z*z
        z4 = z2*z2
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4]) + z4 * (P[5] + z * P[6])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4]) + z4 * (Q[5] + z * Q[6])
        sgn*y*p/q
    end
end
