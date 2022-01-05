"""
    cl5(x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Cl}_5(x)``
for a real angle ``x`` of type `Float64`.  This function is defined as

```math
\\operatorname{Cl}_5(x) = \\Re[\\operatorname{Li}_5(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\cos(kx)}{k^5}
```

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl5(1.0)
0.5228208076420943
```
"""
function cl5(x::Float64)::Float64
    zeta5 = 1.0369277551433699
    pi28 = pi*pi/8.0

    x = range_reduce_odd(x)

    if x == zero(x)
        return zeta5
    end

    if x < 0.5*pi
        P = (1.0369277551433699e+00, -6.1354800479984468e-01,
             9.4076401395712763e-02, -9.4056155866704436e-04)
        Q = (1.0000000000000000e+00, -1.2073698633244778e-02,
             1.3703409625482991e-05, -1.9701280330628469e-09,
             2.1944550184416500e-11)
        y = x*x
        y2 = y*y
        p = P[1] + y * P[2] + y2 * (P[3] + y * P[4])
        q = Q[1] + y * Q[2] + y2 * (Q[3] + y * Q[4] + y2 * Q[5])
        p/q - 1/24*y2*log(x)
    else
        P = (-4.5930112735784898e-01, 4.3720705508867954e-01,
             -7.5895226486465095e-02, 5.2244176912488065e-03,
             -1.5677716622013956e-04, 1.6641624171748576e-06)
        Q = ( 1.0000000000000000e+00, -1.2211486825401188e-01,
              3.8940070749313620e-03, -2.2674805547074318e-05,
             -7.4383354448335299e-08, -3.4131758392216437e-10)
        y = pi - x
        z = y*y - pi28
        z2 = z*z
        z4 = z2*z2
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4]) + z4 * (P[5] + z * P[6])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4]) + z4 * (Q[5] + z * Q[6])
        p/q
    end
end
