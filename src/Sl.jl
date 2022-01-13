# zeta(n) for n = 2,...,9
const zeta = (
    1.6449340668482264, 1.2020569031595943, 1.0823232337111382,
    1.0369277551433699, 1.0173430619844491, 1.0083492773819228,
    1.0040773561979443, 1.0020083928260822
)

"""
    sl(n::Integer, x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Sl}_n(x)``
for integers ``n > 0`` and a real angle ``x`` of type `Float64`.  This
function is defined as

```math
\\begin{aligned}
\\operatorname{Cl}_n(x) &= \\Re[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\cos(kx)}{k^n}, \\qquad \\text{for}~n~\\text{even}, \\\\
\\operatorname{Cl}_n(x) &= \\Im[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\sin(kx)}{k^n}, \\qquad \\text{for}~n~\\text{odd}.
\\end{aligned}
```

The implementation follows the approach presented in [Richard
J. Mathar, "A C99 Implementation of the Clausen Sums", arXiv:1309.7504].

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> sl(10, 1.0)
0.8423605391686301
```
"""
function sl(n::Integer, x::Float64)::Float64
    if n < 1
        throw(DomainError(n, "sl(n,x) undefined for n < 1"))
    end

    (x, sgn) = range_reduce(n + 1, x)

    if n < 10
        coeff = zeros(Float64, n + 1);

        # initialize for n = 1
        coeff[1] = 0.5*pi # coefficient of x^0
        coeff[2] = -0.5   # coefficient of x^1

        # integrate (n - 1) times
        for k in 2:n
            sign_flip = iseven(k) ? -1.0 : 1.0
            # integration
            for i in (k + 1):-1:2
                coeff[i] = sign_flip*coeff[i - 1]/(i - 1)
            end
            # integration constant
            coeff[1] = iseven(k) ? zeta[k - 1] : 0.0
        end

        sgn*evalpoly(x, coeff)
    else
        sgn*sl_series(n, x)
    end
end