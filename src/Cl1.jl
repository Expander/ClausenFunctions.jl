"""
    cl1(x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Cl}_1(x)``
for a real angle ``x`` of type `Float64`.  This function is defined as

```math
\\operatorname{Cl}_1(x) = \\Re[\\operatorname{Li}_1(e^{ix})] = \\Re[-\\log(1 - e^{ix})]
```

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl1(1.0)
0.042019505825369
```
"""
function cl1(x::Float64)::Float64
    if x < zero(x)
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

    if x == zero(x)
        return Inf
    end

    real(-log(one(x) - exp(x*im)))
end
