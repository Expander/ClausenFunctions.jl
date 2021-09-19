is_equal(a::Float64, b::Float64, eps::Float64) =
    abs(a - b) <= (1.0 + max(abs(a), abs(b))) * eps
