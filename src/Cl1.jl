"""
    cl1(x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Cl}_1(x)``
for a real angle ``x`` of type `Float64`.  This function is defined as

```math
\\operatorname{Cl}_1(x) = \\Re[\\operatorname{Li}_1(e^{ix})] = \\Re[-\\log(1 - e^{ix})]
```

Note: ``\\operatorname{Cl}_1(x)`` is not defined for ``x=2n\\pi`` with
``n\\in\\mathbb{Z}``.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl1(1.0)
0.04201950582536895
```
"""
function cl1(x::Float64)::Float64
    x = range_reduce_odd(x)

    if x == zero(x)
        return Inf
    end

    -log(2.0*sin(0.5*x))
end
