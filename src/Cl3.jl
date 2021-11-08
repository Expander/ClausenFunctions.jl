"""
    cl3(x::Float64)::Float64

Returns the value of the Clausen function ``\\mathrm{Cl}_3(x)`` for a
real angle `x` of type `Float64`.  This function is defined as

```math
\mathrm{Cl}_3(x) = \Re \mathrm{Li}_3(e^{ix})
```

Author: Alexander Voigt

License: MIT

# Example
```julia
cl3(1.0)
```
"""
function cl3(x::Float64)::Float64
    zeta3 = 1.2020569031595943
    pi28 = pi*pi/8.0

    if x < 0.0
        x = -x
    end

    if x >= 2.0*pi
        x = mod2pi(x)
    end

    if x > pi
        p0 = 6.28125
        p1 = 0.0019353071795864769253
        x = (p0 - x) + p1
    end

    if x == 0.0
        return zeta3
    end

    if x < 0.5*pi
        P = (-7.5430148591242361e-01,  1.6121940167854339e-02,
             -3.7484056212140535e-05, -2.5191292110169198e-07)
        Q = (1.0000000000000000e+00, -2.6015033560727570e-02,
             1.5460630299236049e-04, -1.0987530650923219e-07)
        y = x*x
        z = y - pi28
        z2 = z*z
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4])
        zeta3 + y*(p/q + 0.5*log(x))
    else
        P = (-4.9017024647634973e-01, 4.1559155224660940e-01,
             -7.9425531417806701e-02, 5.9420152260602943e-03,
             -1.8302227163540190e-04, 1.8027408929418533e-06)
        Q = (1.0000000000000000e+00, -1.9495887541644712e-01,
             1.2059410236484074e-02, -2.5235889467301620e-04,
             1.0199322763377861e-06,  1.9612106499469264e-09)
        y = pi - x
        z = y*y - pi28
        z2 = z*z
        z4 = z2*z2
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4]) + z4 * (P[5] + z * P[6])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4]) + z4 * (Q[5] + z * Q[6])
        p/q
    end
end
