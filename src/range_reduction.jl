# returns 2pi - x, avoiding a catastrophic cancellation
function two_pi_minus(x::Float16)::Float16
    2pi - x
end

# returns 2pi - x, avoiding a catastrophic cancellation
function two_pi_minus(x::Float32)::Float32
    2pi - x
end

# returns 2pi - x, avoiding a catastrophic cancellation
function two_pi_minus(x::Float64)
    p0 = 6.28125
    p1 = 0.0019353071795864769253
    (p0 - x) + p1
end

# Generic fallback for BigFloat and other arbitrary-precision types
function two_pi_minus(x::Real)
    2*oftype(x, pi) - x
end

# returns range-reduced x in [0,pi] for odd n
function range_reduce_odd(x::Real)
    if x < zero(x)
        x = -x
    end

    twopi = 2*oftype(x, pi)

    if x >= twopi
        x = mod(x, twopi)
    end

    if x > oftype(x, pi)
        x = two_pi_minus(x)
    end

    x
end

# returns (x, sign) with range-reduced x in [0,pi] for even n
function range_reduce_even(x::Real)
    sgn = one(x)

    if x < zero(x)
        x = -x
        sgn = -one(x)
    end

    twopi = 2*oftype(x, pi)

    if x >= twopi
        x = mod(x, twopi)
    end

    if x > oftype(x, pi)
        x = two_pi_minus(x)
        sgn = -sgn
    end

    (x, sgn)
end

# returns (x, sign) with range-reduced x in [0,pi]
function range_reduce(n::Integer, x::Real)
    if iseven(n)
        range_reduce_even(x)
    else
        (range_reduce_odd(x), one(x))
    end
end

# returns x - 2πk in [-π, π] for a BigFloat x, rounded to the working
# precision.  The subtraction cancels up to exponent(x) leading bits
# and, near multiples of 2π, up to about precision(x) further bits,
# so 2π is needed with that many extra bits.  MPFR computes the
# remainder exactly before rounding.
function rem_twopi_nearest(x::BigFloat)
    prec = max(precision(x), precision(BigFloat))
    extra = max(exponent(x), 0) + prec + 32
    twopi = setprecision(BigFloat, prec + extra) do
        2*BigFloat(pi)
    end
    rem(x, twopi, RoundNearest)
end

# returns range-reduced x in [0,pi] for odd n
function range_reduce_odd(x::BigFloat)
    isfinite(x) && !iszero(x) || return abs(x - x)
    abs(rem_twopi_nearest(x))
end

# returns (x, sign) with range-reduced x in [0,pi] for even n
function range_reduce_even(x::BigFloat)
    iszero(x) && return (x, one(x))
    isfinite(x) || return (x - x, one(x))
    r = rem_twopi_nearest(x)
    (abs(r), signbit(r) ? -one(x) : one(x))
end

# Branch-free floating-point Payne–Hanek reduction modulo π/2 for
# |x| ≥ 2^14, following T. Ly, "A Performance Improvement of the
# Payne–Hanek Range Reduction Algorithm", arXiv:2609.35015 (Algorithm 2,
# as in LLVM libc), with reduction modulus π/2^N for N = 1, table limbs of
# p_c = 51 bits and one table entry per M = 16 binades.  The reduced
# argument has an absolute error below 2^-110 (relative below 2^-60) for
# every input.
#
# Table entry i holds D₀ + D₁ + D₂ + D₃ ≈ 2^N {2^(M i)/π}, where {⋅} is the
# centered fractional part; inputs x with exponent e use i = ⌊(e - 62)/M⌋
# and are scaled exactly to x_r = x 2^(-M i).  Entries i < 0 carry the
# prefix 2^(N-2) instead, whose product with x_r vanishes modulo 2^(N+1).
const _PH_N = 1
const _PH_M = 16
const _PH_IMIN = -3   # covers |x| ≥ 2^14
const _PH_IMAX = 60   # covers |x| < 2^1024

function _ph_table()
    setprecision(BigFloat, 2048) do
        limb(v, bits) = Float64(BigFloat(v; precision = bits))
        map(_PH_IMIN:_PH_IMAX) do i
            a = big(2.0)^(_PH_M*i)/BigFloat(pi)
            v = big(2.0)^_PH_N*(i >= 0 ? a - round(a) : big(0.25) + a)
            D0 = limb(v, 51); v -= D0
            D1 = limb(v, 51); v -= D1
            D2 = limb(v, 51); v -= D2
            (D0, D1, D2, Float64(v))
        end
    end
end

const _PH_TABLE = _ph_table()

# returns (n mod 4, yh, yl) with x = n π/2 + (yh + yl), |yh + yl| ≲ π/4,
# for finite x with |x| ≥ 2^14
@inline function rem_pio2_large(x::Float64)
    i = (exponent(x) - 62) >> 4                              # ⌊(e - 62)/16⌋
    xr = x*reinterpret(Float64, UInt64(1023 - _PH_M*i) << 52)  # exact
    @inbounds D0, D1, D2, D3 = _PH_TABLE[i - _PH_IMIN + 1]
    # exact products; x_r D₀'s high part is ≡ 0 modulo 4 and is dropped
    ah = xr*D0; al = fma(xr, D0, -ah)
    bh = xr*D1; bl = fma(xr, D1, -bh)
    ch = xr*D2; cl = fma(xr, D2, -ch)
    k = round(al + bh)          # quotient
    v = (al - k) + bh           # exact
    qh = bl + ch                # Fast2Sum, exact
    ql = ch - (qh - bl)
    r = fma(xr, D3, cl)
    uh = v + qh                 # Fast2Sum, exact
    ul = (qh - (uh - v)) + (ql + r)
    # u = x 2/π - k as double-double; y = u π/2, renormalized (Fast2Sum)
    # so that |yl| ≤ ulp(yh)/2
    ph = uh*1.5707963267948966
    pl = fma(uh, 1.5707963267948966, -ph) + fma(uh, 6.123233995736766e-17, ul*1.5707963267948966)
    yh = ph + pl
    yl = pl - (yh - ph)
    (unsafe_trunc(Int, k) & 3, yh, yl)
end
