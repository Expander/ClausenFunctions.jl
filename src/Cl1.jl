"""
    cl1(x::Real)::Real

Returns the value of the Clausen function ``\\operatorname{Cl}_1(x)``
for a real angle ``x`` of type `Real`.  This function is defined as

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
function cl1(x::Real)::Real
    x = range_reduce_odd(x)

    if x == zero(x)
        return Inf
    end

    -log(2*sin(one(x)/2*x))
end
