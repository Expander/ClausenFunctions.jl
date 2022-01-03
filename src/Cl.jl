import SpecialFunctions

# returns (x, sign) with x in [0,pi]
function range_reduce(n::UInt64, x::Float64)
    sgn = one(Float64)

    if x < zero(Float64)
        x = -x
        sgn = -one(Float64)
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

    if iseven(n)
        (x, sgn)
    else
        (x, one(Float64))
    end
end

# returns B_{2k}/(2k)! where B_{2k} are the even Bernoulli numbers
function b2(k::UInt64)
    2*(-1)^(k + 1)*SpecialFunctions.zeta(2*k)/(2*pi)^(2*k)
end

# Eq.(2.11)
function ncal(n::UInt64, x::Float64)
    sum = zero(Float64)
    k = one(UInt64)

    while k < typemax(UInt64)
        term = (-1)^k*b2(k)*x^(2*k + n + 1)/((2*k + n + 1))
        old_sum = sum
        sum += term
        sum == old_sum && break
        k += 1
    end

    (x^(n + 1)/(n + 1) + sum)/(n + 1)
end

# returns Cl(n, 0)
function cl0(n::UInt64)
    if iseven(n)
        zero(Float64)
    else
        SpecialFunctions.zeta(n)
    end
end

function pcal(k::UInt64, x::Float64)
    sum = zero(Float64)

    for i in 2:k
       sum += (-1)^(fld(k - 1, 2.0) + fld(i - 1, 2.0))*x^(k - i)/factorial(k - i)*cl0(i)
    end

    sum
end

"""
    cl(n::UInt64, x::Float64)::Float64

Returns the value of the Clausen function ``\\operatorname{Cl}_n(x)``
for integer `n > 1` and a real angle `x` of type `Float64`.  This
function is defined as

```math
\\operatorname{Cl}_n(x) &= \\Im[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\sin(kx)}{k^n}, \\qquad n~\\text{even}, \\
\\operatorname{Cl}_n(x) &= \\Re[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\cos(kx)}{k^n}, \\qquad n~\\text{odd}
```

The implementation follows the approach presented in [Jiming Wu,
Xiaoping Zhang, Dongjie Liu, "An efficient calculation of the Clausen
functions Cl_n(0)(n >= 2)"].

Author: Alexander Voigt

License: MIT

# Example
```julia
cl(10, 1.0)
```
"""
function cl(n::UInt64, x::Float64)::Float64
    if n < 2
        throw("cl(n,x) undefined for n < 2")
    end

    (x, sgn) = range_reduce(n, x)

    if x == zero(Float64)
        return cl0(n)
    end

    if iseven(n) && x == pi
        return zero(Float64)
    end

    # sum in Eq.(2.8), (2.13)
    sum = zero(Float64)

    for i in 0:(n - 2)
       sum += one(Float64)*(-1)^i*binomial(n - 2, i)*x^i*ncal(n - 2 - i, x)
    end

    # Eq.(2.8), (2.13)
    sgn*((-1)^fld(n + 1, 2.0)*x^(n - 1)/factorial(n - 1)*log(2*sin(x/2))
         + (-1)^(fld(n, 2.0) + 1)/factorial(n - 2)*sum + pcal(n, x))
end
