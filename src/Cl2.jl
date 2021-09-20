"""
    cl2(x::Float64)::Float64

Returns the value of the Clausen function ``\\mathrm{Cl}_2(x)`` for a
real angle `x` of type `Float64`.  The Clausen function is defined as
```\\mathrm{Cl}_2(x) = -\\int_0^x \\log|2\\sin(t/2)| dt```

Author: Alexander Voigt

License: MIT

# Example
```julia
cl2(1.0)
```
"""
function cl2(x::Float64)::Float64
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
        P = (
            2.7951565822419270e-02, -8.8865360514541522e-04,
            6.8282348222485902e-06, -7.5276232403566808e-09
        )
        Q = (
            1.0000000000000000e+00, -3.6904397961160525e-02,
            3.7342870576106476e-04, -8.7460760866531179e-07
        )
        y = x*x
        z = y - pi28
        z2 = z*z
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4])
        sgn*x*(1.0 - log(x) + y*p/q/2.0)
    else
        P = (
            6.4005702446195512e-01, -2.0641655351338783e-01,
            2.4175305223497718e-02, -1.2355955287855728e-03,
            2.5649833551291124e-05, -1.4783829128773320e-07
        )
        Q = (
            1.0000000000000000e+00, -2.5299102015666356e-01,
            2.2148751048467057e-02, -7.8183920462457496e-04,
            9.5432542196310670e-06, -1.8184302880448247e-08
        )
        y = pi - x
        z = y*y - pi28
        z2 = z*z
        z4 = z2*z2
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4]) + z4 * (P[5] + z * P[6])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4]) + z4 * (Q[5] + z * Q[6])
        sgn*y*p/q
    end
end
