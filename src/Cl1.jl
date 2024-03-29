"""
    cl1(z)

Returns the value of the Clausen function ``\\operatorname{Cl}_1(z)``
for an argument ``z`` of type `Real` or `Complex`.  This function is
defined as

```math
\\operatorname{Cl}_1(z) = -\\frac{1}{2}\\left[\\log(1 - e^{iz}) + \\log(1 - e^{-iz})\\right]
```

Note: ``\\operatorname{Cl}_1(z)`` is not defined for ``z=2k\\pi`` with
``k\\in\\mathbb{Z}``.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl1(1.0)
0.04201950582536895

julia> cl1(big"1")
0.04201950582536896172579838403790203712453892055703441769956888996856898991572203

julia> cl1(1.0 + 1.0im)
-0.3479608285425303 - 0.7021088550913619im

julia> cl1(big"1" + 1im)
-0.3479608285425304681676697311626225254052481291939705626677618712960309352197459 - 0.7021088550913618982002640571209884840333848751442483836897991093694470903764252im
```
"""
cl1

cl1(x::Float16) = oftype(x, cl1(Float32(x)))

cl1(x::Float32) = oftype(x, cl1(Float64(x)))

function cl1(x::Float64)::Float64
    x = range_reduce_odd(x)

    if iszero(x)
        return Inf
    end

    -log(2*sin(one(x)/2*x))
end

function cl1(x::Real)
    -log(hypot(one(x) - cos(x), sin(x)))
end

function cl1(z::Complex)
    eiz = cis(z)
    -(log(one(z) - eiz) + log(one(z) - inv(eiz)))/2
end
