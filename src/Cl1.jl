"""
    cl1(x::Real)::Real

Returns the value of the Clausen function ``\\operatorname{Cl}_1(x)``
for a real angle ``x`` of type `Real`.  This function is defined as

```math
\\operatorname{Cl}_1(x) = \\Re[\\operatorname{Li}_1(e^{ix})] = \\Re[-\\log(1 - e^{ix})]
```

Note: ``\\operatorname{Cl}_1(x)`` is not defined for ``x=2k\\pi`` with
``k\\in\\mathbb{Z}``.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl1(1.0)
0.04201950582536895
```
"""
cl1(x::Real) = _cl1(float(x))

_cl1(x::Float16) = oftype(x, _cl1(Float32(x)))

_cl1(x::Float32) = oftype(x, _cl1(Float64(x)))

function _cl1(x::Float64)::Float64
    x = range_reduce_odd(x)

    if x == zero(x)
        return Inf
    end

    -log(2*sin(one(x)/2*x))
end
