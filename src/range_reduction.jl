# returns range-reduced x in [0,pi] for odd n
function range_reduce_odd(x::Real)
    if x < zero(x)
        x = -x
    end

    if x >= 2pi
        x = mod(x, 2pi)
    end

    if x > pi
        p0 = 6.28125
        p1 = 0.0019353071795864769253
        x = (p0 - x) + p1
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

    if x >= 2pi
        x = mod(x, 2pi)
    end

    if x > pi
        p0 = 6.28125
        p1 = 0.0019353071795864769253
        x = (p0 - x) + p1
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
