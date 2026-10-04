# Generates the minimax polynomials _CL2_P1, _CL2_P2 and _CL2_P2_LO
# used by cl2(::Float64) in src/Cl2.jl.
#
# Usage (Remez.jl is not a dependency of ClausenFunctions):
#
#   julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("Remez"); include("math/cl2_minimax.jl")'

using Remez

setprecision(BigFloat, 512)

const SPLIT = big(pi)/3 # branch point between the two expansions
const DEG1 = 6          # degree of P₁
const DEG2 = 16         # degree of P₂
const NTERMS = 400      # Taylor terms used to evaluate the targets

# Bernoulli numbers B_0, ..., B_N as exact rationals (Akiyama-Tanigawa)
function bernoulli(N)
    A = Vector{Rational{BigInt}}(undef, N + 1)
    B = similar(A)
    for m in 0:N
        A[m + 1] = 1//(m + 1)
        for j in m:-1:1
            A[j] = j*(A[j] - A[j + 1])
        end
        B[m + 1] = A[1]
    end
    B # B[n + 1] = B_n
end

const BN = bernoulli(2NTERMS)

horner(t, c) = foldr((ci, s) -> muladd(s, t, ci), c; init = zero(t))

# Cl₂(a) = a - a log(a) + a³ g₁(a²) with
# g₁(t) = Σ_{n≥1} |B_{2n}|/(2n (2n+1) (2n)!) t^(n-1)
const C1 = [BigFloat(abs(BN[2n + 1]))/(2n*(2n + 1)*factorial(big(2n))) for n in 1:NTERMS]
g1(t) = horner(t, C1)
cl2_near0(a) = a - a*log(a) + a^3*g1(a^2)

# Cl₂(π - y) = y h(y²) with h(t) = log(2) + t g₂(t),
# g₂(t) = Σ_{m≥1} (-1)^m (2^(2m) - 1) B_{2m}/(2m (2m+1)!) t^(m-1)
const C2 = [BigFloat((-1)^m*(big(2)^(2m) - 1)*BN[2m + 1])/(2m*factorial(big(2m + 1))) for m in 1:NTERMS]
g2(t) = horner(t, C2)
h(t) = log(big(2)) + t*g2(t)

# weights: error relative to Cl₂ itself
w1(t, _) = (a = sqrt(t); iszero(t) ? zero(t) : a^3/cl2_near0(a))
w2(t, _) = t/h(t)

P1, _, E1, _ = ratfn_minimax(g1, (big(0), SPLIT^2), DEG1, 0, w1)
P2, _, E2, _ = ratfn_minimax(g2, (big(0), (big(pi) - SPLIT)^2), DEG2, 0, w2)

println("# relative errors: P₁ ", Float64(E1), ", P₂ ", Float64(E2))
println("const _CL2_P1 = ", Tuple(Float64.(P1)))
println("const _CL2_P2_LO = ", Float64(P2[1] - Float64(P2[1])))
println("const _CL2_P2 = ", Tuple(Float64.(P2)))
