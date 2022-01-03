"""
    cl6(x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Cl}_n(x)``
for positive integer `n` and a real angle `x` of type `Float64`.  This
function is defined as

```math
\\operatorname{Cl}_n(x) &= \\Im[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\sin(kx)}{k^n}, \\qquad n~\\text{even}, \\
\\operatorname{Cl}_n(x) &= \\Re[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\cos(kx)}{k^n}, \\qquad n~\\text{odd}
```

Author: Alexander Voigt

License: MIT

# Example
```julia
cl(10, 1.0)
```
"""
function cl(n::UInt64, x::Float64)::Float64
    0.0
end
