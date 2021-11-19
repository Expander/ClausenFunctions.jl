"""
    cl4(x::Float64)::Float64

Returns the value of the Clausen function ``\\mathrm{Cl}_4(x)`` for a
real angle `x` of type `Float64`.  This function is defined as

```math
\\mathrm{Cl}_4(x) = \\Im[\\mathrm{Li}_4(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\sin(kx)}{k^4}
```

Author: Alexander Voigt

License: MIT

# Example
```julia
cl4(1.0)
```
"""
function cl4(x::Float64)::Float64
    zeta3 = 1.2020569031595943
    pi28 = pi*pi/8.0
    sgn = 1.0

    if x < 0.0
        x = -x
        sgn = -1.0
    end

    if x >= 2.0*pi
        x = mod2pi(x)
    end

    if x > pi
        p0 = 6.28125
        p1 = 0.0019353071795864769253
        x = (p0 - x) + p1
        sgn = -sgn
    end

    if x == 0.0 || x == pi
        return 0.0
    end

    if x < 0.5*pi
        P = (-3.0641482939025622e-01,  6.1700119337868541e-03,
             -2.0244370294666391e-05, -3.1997168417491939e-08)
        Q = (1.0000000000000000e+00, -2.2415973860234228e-02,
             1.1164184654978598e-04, -6.3541742491831717e-08)
        y = x*x
        z = y - pi28
        z2 = z*z
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4])
        sgn*x*(zeta3 + y*(p/q + log(x)/6))
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
