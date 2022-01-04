"""
    cl3(x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Cl}_3(x)``
for a real angle `x` of type `Float64`.  This function is defined as

```math
\\operatorname{Cl}_3(x) = \\Re[\\operatorname{Li}_3(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\cos(kx)}{k^3}
```

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl3(1.0)
0.44857300728001737
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
        P = (-7.5000000000000001e-01,  1.5707637881835541e-02,
             -3.5426736843494423e-05, -2.4408931585123682e-07)
        Q = (1.0000000000000000e+00, -2.5573146805410089e-02,
             1.5019774853075050e-04, -1.0648552418111624e-07)
        y = x*x
        y2 = y*y
        p = P[1] + y * P[2] + y2 * (P[3] + y * P[4])
        q = Q[1] + y * Q[2] + y2 * (Q[3] + y * Q[4])
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
