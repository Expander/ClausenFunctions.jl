"""
    cl2(x::Real)::Real

Returns the value of the Clausen function ``\\operatorname{Cl}_2(x)``
for a real angle ``x`` of type `Real`.  The Clausen function is
defined as

```math
\\operatorname{Cl}_2(x) = -\\int_0^x \\log|2\\sin(t/2)| dt
```

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> cl2(1.0)
1.0139591323607684
```
"""
cl2(x::Real) = _cl2(float(x))

_cl2(x::Float16) = oftype(x, _cl2(Float32(x)))

_cl2(x::Float32) = oftype(x, _cl2(Float64(x)))

function _cl2(x::Float64)::Float64
    (x, sgn) = range_reduce_even(x)

    if iszero(x)
        x
    elseif x < pi/2
        P = (1.3888888888888889e-02, -4.3286930203743071e-04,
             3.2779814789973427e-06, -3.6001540369575084e-09)
        Q = (1.0000000000000000e+00, -3.6166589746694121e-02,
             3.6015827281202639e-04, -8.3646182842184428e-07)
        y = x*x
        y2 = y*y
        p = P[1] + y * P[2] + y2 * (P[3] + y * P[4])
        q = Q[1] + y * Q[2] + y2 * (Q[3] + y * Q[4])
        sgn*x*(1.0 - log(x) + y*p/q)
    elseif x < pi
        P = (6.4005702446195512e-01, -2.0641655351338783e-01,
             2.4175305223497718e-02, -1.2355955287855728e-03,
             2.5649833551291124e-05, -1.4783829128773320e-07)
        Q = (1.0000000000000000e+00, -2.5299102015666356e-01,
             2.2148751048467057e-02, -7.8183920462457496e-04,
             9.5432542196310670e-06, -1.8184302880448247e-08)
        y = pi - x
        z = y*y - pi*pi/8
        z2 = z*z
        z4 = z2*z2
        p = P[1] + z * P[2] + z2 * (P[3] + z * P[4]) + z4 * (P[5] + z * P[6])
        q = Q[1] + z * Q[2] + z2 * (Q[3] + z * Q[4]) + z4 * (Q[5] + z * Q[6])
        sgn*y*p/q
    end
end


# BigFloat Cl₂ via Kummer-accelerated Bernoulli series:
#
# Cl₂(θ)/θ = 3 - log[θ(1 - θ²/(4π²))] - (2π/θ)log((2π+θ)/(2π-θ))
#           + Σ_{n=1}^K (ζ(2n)-1)/(n(2n+1))  (θ/(2π))^{2n}
#
# valid for |θ| < 2π (always true after range reduction to [0, π]).
#
# Rewriting in the θ^{2n} basis avoids large intermediate values:
#
#   d_n = (ζ(2n)-1)/(n(2n+1)(2π)^{2n})
#       = |B_{2n}|/(2n(2n+1)(2n)!)  -  1/(n(2n+1)(2π)^{2n})
#
# The first term is exactly the standard Bernoulli coefficient (cheap,
# small values); the second is a trivial geometric correction.
#
# Compared to the unaccelerated series, ζ(2n)-1 ≈ 1/4^n roughly halves
# K at the worst case θ = π, while the coefficient computation stays
# as cheap as the original.
#
# Convergence rate: (θ/(4π))² per term (for large n).
# At θ = π (worst case): ~1.2 digits per term → K ≈ 0.83P for P digits.
# At θ = π/4: ~2.4 digits per term → K ≈ 0.42P.

# Cache: stores d_n for n = 1, ..., K.
const _cl2_coeff_cache = Ref{Vector{BigFloat}}(BigFloat[])
const _cl2_coeff_prec = Ref{Int}(0)
const _cl2_coeff_K = Ref{Int}(0)
const _cl2_coeff_lock = ReentrantLock()


function _cl2(x::BigFloat)
    (x, sgn) = range_reduce_even(x)

    iszero(x) && return x
    x == BigFloat(pi) && return zero(x)

    # +10: guard digits to absorb rounding in the Horner loop and prefix.
    target_digits = ceil(Int, precision(BigFloat)*log10(2)) + 10

    twopi = 2*BigFloat(pi)
    v = x*x

    # Each series term decays by a factor of u = (θ/(2π))², giving
    # -2·log₁₀(θ/(2π)) decimal digits per term.  The Kummer acceleration
    # replaces ζ(2n) with ζ(2n)-1 ≈ 1/2^{2n} = 1/4ⁿ, contributing an
    # extra log₁₀(4) ≈ 0.602 digits per term (rounded down to 0.6).
    # The floor of 0.5 prevents division by zero for θ near 2π (which
    # cannot occur after range reduction to [0, π], but is defensive).
    ratio = Float64(x)/(2*Float64(pi))
    digits_per_term = max(-2*log10(ratio) + 0.6, 0.5)
    # +30: safety margin covering (a) the first ~10 terms where
    # ζ(2n)-1 ≫ 1/4ⁿ so convergence is slower than the asymptotic
    # estimate, and (b) any rounding in the Float64 K estimate.
    K = ceil(Int, target_digits/digits_per_term) + 30

    c = _cl2_ensure_coeffs(K)

    # Horner: S = v(d₁ + v(d₂ + ... + v dₖ))
    s = c[K]
    for k in (K - 1):-1:1
         s = muladd(s, v, c[k])
    end
    s *= v

    # Prefix: 3 - log[θ(1 - θ²/(4π²))] - (2π/θ)log((2π+θ)/(2π-θ))
    # log1p avoids cancellation for small θ where (2π+θ)/(2π-θ) ≈ 1.
    prefix = 3 - log(x*(one(x) - v/twopi^2)) - (twopi/x)*log1p(2*x/(twopi - x))

    sgn*x*(prefix + s)
end


# Compute Bernoulli numbers B_0, B_2, ..., B_{2K} via the standard recurrence.
function _bernoulli_even(K::Int)
    n = 2K
    B = zeros(BigFloat, n + 1)
    B[1] = one(BigFloat)
    B[2] = -one(BigFloat)/2

    for m in 2:n
        s = zero(BigFloat)
        binom = one(BigFloat)
        for k in 0:m-1
            s = binom*B[k + 1] + s
            binom = binom*(m + 1 - k)/(k + 1)
        end
        B[m + 1] = -s/(m + 1)
    end

    Beven = Vector{BigFloat}(undef, K + 1)

    for i in 0:K
        Beven[i + 1] = B[2i + 1]
    end

    Beven
end


# Compute and cache d_n = |B_{2n}|/(2n(2n+1)(2n)!) - 1/(n(2n+1)(2π)^{2n})
# for n = 1,...,K.
function _cl2_ensure_coeffs(K::Int)
    prec = precision(BigFloat)

    lock(_cl2_coeff_lock) do
        if _cl2_coeff_prec[] >= prec && _cl2_coeff_K[] >= K
            _cl2_coeff_cache[]
        else
            Beven = _bernoulli_even(K)

            inv_twopi_sq = 1/(2*BigFloat(pi))^2
            c = Vector{BigFloat}(undef, K)
            fac = one(BigFloat)
            inv_twopi_pow = one(BigFloat)

            for k in 1:K
                m = 2k
                fac *= (m - 1)*m # (2k)!
                inv_twopi_pow *= inv_twopi_sq # 1/(2π)^{2k}
                # Standard Bernoulli coefficient: |B_{2k}|/(2k(2k+1)(2k)!)
                orig = abs(Beven[k + 1])/(m*(m + 1)*fac)
                # Kummer correction: 1/(k(2k+1)(2π)^{2k})
                corr = inv_twopi_pow/(k*(m + 1))
                c[k] = orig - corr
            end

            _cl2_coeff_cache[] = c
            _cl2_coeff_prec[] = prec
            _cl2_coeff_K[] = K

            c
        end
    end
end
