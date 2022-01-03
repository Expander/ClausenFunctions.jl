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

# bernoulli numbers
function bernoulli2(n::UInt64)
    if n > 35
        throw(DomainError(n, "If n > 35, then the numerator needs Int128 at least, and worse, so this code is not the code you want. Try using bernoulli(n, 0.0) to get a floating point approximation to the result."))
    end

    # Denominator of Bernoulli number B_n
    #   http://oeis.org/A027642
    D = [2, 6, 1, 30, 1, 42, 1, 30, 1, 66, 1, 2730, 1, 6, 1, 510, 1, 798, 1, 330, 1, 138, 1, 2730, 1, 6, 1, 870, 1, 14322, 1, 510, 1, 6, 1, 1919190, 1, 6, 1, 13530, 1, 1806, 1, 690, 1, 282, 1, 46410, 1, 66, 1, 1590, 1, 798, 1, 870, 1, 354, 1, 56786730]

    # Numerator of Bernoulli number B_n (storing 62 of these because they are easy)
    #   http://oeis.org/A027641
    N = [-1, 1, 0, -1, 0, 1, 0, -1, 0, 5, 0, -691, 0, 7, 0, -3617, 0, 43867, 0, -174611, 0, 854513, 0, -236364091, 0, 8553103, 0, -23749461029, 0, 8615841276005, 0, -7709321041217, 0, 2577687858367, 1]

    N[n] // D[n]
end

# Eq.(2.11)
function ncal(n::UInt64, x::Float64)
    sum = zero(Float64)
    k = one(UInt64)

    while k < 17
        term = (-1)^k*bernoulli2(2*k)*x^(2*k + n + 1)/((2*k + n + 1)*factorial(big(2*k)))
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
