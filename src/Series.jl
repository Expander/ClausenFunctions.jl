# returns sum(k=1:inf, cos(k*x)/k^n)
function co_series(n::Integer, x::Number)
    kmax = max_terms(n, typeof(x))
    co = cos(x)
    co2 = one(x) # cos((k-2)*x)
    co1 = co     # cos((k-1)*x)
    sum = co
    for k in 2:kmax
        con = 2*co*co1 - co2 # cos(k*x)
        co2 = co1
        co1 = con
        sum += con/oftype(x, k)^n
    end
    sum
end

# returns sum(k=1:inf, sin(k*x)/k^n)
function si_series(n::Integer, x::Number)
    kmax = max_terms(n, typeof(x))
    si, co = sincos(x)
    si2 = zero(x) # sin((k-2)*x)
    si1 = si      # sin((k-1)*x)
    sum = si
    for k in 2:kmax
        si = 2*co*si1 - si2 # sin(k*x)
        si2 = si1
        si1 = si
        sum += si/oftype(x, k)^n
    end
    sum
end

# returns maximum number of terms in co_series or si_series
function max_terms(n::Integer, ::Type{T}) where T
    ceil(typeof(n), eps(T)^(-inv(n)))
end
