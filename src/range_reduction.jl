# returns (x, sign) with x in [0,pi]
function range_reduce(n::Integer, x::Float64)
    sgn = one(x)

    if x < zero(x)
        x = -x
        sgn = -one(x)
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
        (x, one(x))
    end
end
