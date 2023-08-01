# returns Sl(n,x) using the naive series expansion
function sl_series(n::Integer, x::Number)
    kmax = ceil(typeof(n), eps(typeof(x))^(-inv(n)))

    if iseven(n)
        co = cos(x)
        co2 = one(x) # cos((n-2)*x)
        co1 = co     # cos((n-1)*x)
        sum = co
        for k in 2:kmax
            con = 2*co*co1 - co2 # cos(n*x)
            co2 = co1
            co1 = con
            sum += con/oftype(x, k)^n
        end
        sum
    else
        si, co = sincos(x)
        si2 = zero(x) # sin((n-2)*x)
        si1 = si      # sin((n-1)*x)
        sum = si
        for k in 2:kmax
            si = 2*co*si1 - si2 # sin(n*x)
            si2 = si1
            si1 = si
            sum += si/oftype(x, k)^n
        end
        sum
    end
end

"""
    sl(n::Integer, x::Real)::Real

Returns the value of the Glaisher-Clausen function
``\\operatorname{Sl}_n(x)`` for integers ``n > 0`` and a real angle
``x`` of type `Real`.  This function is defined as

```math
\\begin{aligned}
\\operatorname{Sl}_n(x) &= \\Re[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\cos(kx)}{k^n}, \\qquad \\text{for}~n~\\text{even}, \\\\
\\operatorname{Sl}_n(x) &= \\Im[\\operatorname{Li}_n(e^{ix})] = \\sum_{k=1}^\\infty \\frac{\\sin(kx)}{k^n}, \\qquad \\text{for}~n~\\text{odd}.
\\end{aligned}
```

Note: We set ``\\operatorname{Sl}_1(0) = 0`` for consistency with the
series expansion.

Author: Alexander Voigt

License: MIT

# Example
```jldoctest; setup = :(using ClausenFunctions), output = false
julia> sl(10, 1.0)
0.5398785706335891
```
"""
sl(n::Integer, x::Real) = _sl(n, float(x))

_sl(n::Integer, x::Float16) = oftype(x, _sl(n, Float32(x)))

_sl(n::Integer, x::Float32) = oftype(x, _sl(n, Float64(x)))

function _sl(n::Integer, x::Float64)::Float64
    n < 1 && throw(DomainError(n, "sl(n,x) undefined for n < 1"))

    (x, sgn) = range_reduce(n + 1, x)

    n ==  1 && x == zero(x) && return zero(x)
    n ==  1 && return sgn*(pi/2 - 1/2*x)
    n ==  2 && return sgn*(pi^2/6 + (-1/2*pi + 1/4*x)*x)
    n ==  3 && return sgn*(x*(pi^2/6 + (-1/4*pi + 1/12*x)*x))

    x2 = x*x

    n ==  4 && return sgn*(pi^4/90 + x2*(-1/12*pi^2 + (pi/12 - 1/48*x)*x))
    n ==  5 && return sgn*(x*(pi^4/90 + x2*(-1/36*pi^2 + (pi/48 - 1/240*x)*x)))
    n ==  6 && return sgn*(pi^6/945 + x2*(-1/180*pi^4 + x2*(pi^2/144 + (-1/240*pi + 1/1440*x)*x)))
    n ==  7 && return sgn*(x*(pi^6/945 + x2*(-1/540*pi^4 + x2*(pi^2/720 + (-1/1440*pi + 1/10080*x)*x))))
    n ==  8 && return sgn*(pi^8/9450 + x2*(-1/1890*pi^6 + x2*(pi^4/2160 + x2*(-1/4320*pi^2 + (pi/10080 - 1/80640*x)*x))))
    n ==  9 && return sgn*(x*(pi^8/9450 + x2*(-1/5670*pi^6 + x2*(pi^4/10800 + x2*(-1/30240*pi^2 + (pi/80640 - 1/725760*x)*x)))))
    n == 10 && return sgn*(pi^10/93555 + x2*(-1/18900*pi^8 + x2*(pi^6/22680 + x2*(-1/64800*pi^4 + x2*(pi^2/241920 + (-1/725760*pi + 1/7257600*x)*x)))))
    n == 11 && return sgn*(x*(pi^10/93555 + x2*(-1/56700*pi^8 + x2*(pi^6/113400 + x2*(-1/453600*pi^4 + x2*(pi^2/2177280 + (-1/7257600*pi + 1/79833600*x)*x))))))
    n == 12 && return sgn*((691*pi^12)/638512875 + x2*(-1/187110*pi^10 + x2*(pi^8/226800 + x2*(-1/680400*pi^6 + x2*(pi^4/3628800 + x2*(-1/21772800*pi^2 + (pi/79833600 - 1/958003200*x)*x))))))
    n == 13 && return sgn*(x*((691*pi^12)/638512875 + x2*(-1/561330*pi^10 + x2*(pi^8/1134000 + x2*(-1/4762800*pi^6 + x2*(pi^4/32659200 + x2*(-1/239500800*pi^2 + (pi/958003200 - 1/12454041600*x)*x)))))))
    n == 14 && return sgn*((2*pi^14)/18243225 + x2*((-691*pi^12)/1277025750 + x2*(pi^10/2245320 + x2*(-1/6804000*pi^8 + x2*(pi^6/38102400 + x2*(-1/326592000*pi^4 + x2*(pi^2/2874009600 + (-1/12454041600*pi + 1/174356582400*x)*x)))))))
    n == 15 && return sgn*(x*((2*pi^14)/18243225 + x2*((-691*pi^12)/3831077250 + x2*(pi^10/11226600 + x2*(-1/47628000*pi^8 + x2*(pi^6/342921600 + x2*(-1/3592512000*pi^4 + x2*(pi^2/37362124800 + (-1/174356582400*pi + 1/2615348736000*x)*x))))))))
    n == 16 && return sgn*((3617*pi^16)/325641566250 + x2*(-1/18243225*pi^14 + x2*((691*pi^12)/15324309000 + x2*(-1/67359600*pi^10 + x2*(pi^8/381024000 + x2*(-1/3429216000*pi^6 + x2*(pi^4/43110144000 + x2*(-1/523069747200*pi^2 + (pi/2615348736000 - 1/41845579776000*x)*x))))))))
    n == 17 && return sgn*(x*((3617*pi^16)/325641566250 + x2*(-1/54729675*pi^14 + x2*((691*pi^12)/76621545000 + x2*(-1/471517200*pi^10 + x2*(pi^8/3429216000 + x2*(-1/37721376000*pi^6 + x2*(pi^4/560431872000 + x2*(-1/7846046208000*pi^2 + (pi/41845579776000 - 1/711374856192000*x)*x)))))))))
    n == 18 && return sgn*((43867*pi^18)/38979295480125 + x2*((-3617*pi^16)/651283132500 + x2*(pi^14/218918700 + x2*((-691*pi^12)/459729270000 + x2*(pi^10/3772137600 + x2*(-1/34292160000*pi^8 + x2*(pi^6/452656512000 + x2*(-1/7846046208000*pi^4 + x2*(pi^2/125536739328000 + (-1/711374856192000*pi + 1/12804747411456000*x)*x)))))))))
    n == 19 && return sgn*(x*((43867*pi^18)/38979295480125 + x2*((-3617*pi^16)/1953849397500 + x2*(pi^14/1094593500 + x2*((-691*pi^12)/3218104890000 + x2*(pi^10/33949238400 + x2*(-1/377213760000*pi^8 + x2*(pi^6/5884534656000 + x2*(-1/117690693120000*pi^4 + x2*(pi^2/2134124568576000 + (-1/12804747411456000*pi + 1/243290200817664000*x)*x))))))))))
    n == 20 && return sgn*((174611*pi^20)/1531329465290625 + x2*((-43867*pi^18)/77958590960250 + x2*((3617*pi^16)/7815397590000 + x2*(-1/6567561000*pi^14 + x2*((691*pi^12)/25744839120000 + x2*(-1/339492384000*pi^10 + x2*(pi^8/4526565120000 + x2*(-1/82383485184000*pi^6 + x2*(pi^4/1883051089920000 + x2*(-1/38414242234368000*pi^2 + (pi/243290200817664000 - 1/4865804016353280000*x)*x))))))))))
    n == 21 && return sgn*(x*((174611*pi^20)/1531329465290625 + x2*((-43867*pi^18)/233875772880750 + x2*((3617*pi^16)/39076987950000 + x2*(-1/45972927000*pi^14 + x2*((691*pi^12)/231703552080000 + x2*(-1/3734416224000*pi^10 + x2*(pi^8/58845346560000 + x2*(-1/1235752277760000*pi^6 + x2*(pi^4/32011868528640000 + x2*(-1/729870602452992000*pi^2 + (pi/4865804016353280000 - 1/102181884343418880000*x)*x)))))))))))
    n == 22 && return sgn*((155366*pi^22)/13447856940643125 + x2*((-174611*pi^20)/3062658930581250 + x2*((43867*pi^18)/935503091523000 + x2*((-3617*pi^16)/234461927700000 + x2*(pi^14/367783416000 + x2*((-691*pi^12)/2317035520800000 + x2*(pi^10/44812994688000 + x2*(-1/823834851840000*pi^8 + x2*(pi^6/19772036444160000 + x2*(-1/576213633515520000*pi^4 + x2*(pi^2/14597412049059840000 + (-1/102181884343418880000*pi + 1/2248001455555215360000*x)*x)))))))))))
    n == 23 && return sgn*(x*((155366*pi^22)/13447856940643125 + x2*((-174611*pi^20)/9187976791743750 + x2*((43867*pi^18)/4677515457615000 + x2*((-3617*pi^16)/1641233493900000 + x2*(pi^14/3310050744000 + x2*((-691*pi^12)/25487390728800000 + x2*(pi^10/582568930944000 + x2*(-1/12357522777600000*pi^8 + x2*(pi^6/336124619550720000 + x2*(-1/10948059036794880000*pi^4 + x2*(pi^2/306545653030256640000 + (-1/2248001455555215360000*pi + 1/51704033477769953280000*x)*x))))))))))))
    n == 24 && return sgn*((236364091*pi^24)/201919571963756521875 + x2*((-77683*pi^22)/13447856940643125 + x2*((174611*pi^20)/36751907166975000 + x2*((-43867*pi^18)/28065092745690000 + x2*((3617*pi^16)/13129867951200000 + x2*(-1/33100507440000*pi^14 + x2*((691*pi^12)/305848688745600000 + x2*(-1/8155965033216000*pi^10 + x2*(pi^8/197720364441600000 + x2*(-1/6050243151912960000*pi^6 + x2*(pi^4/218961180735897600000 + x2*(-1/6744004366665646080000*pi^2 + (pi/51704033477769953280000 - 1/1240896803466478878720000*x)*x))))))))))))
    n == 25 && return sgn*(x*((236364091*pi^24)/201919571963756521875 + x2*((-77683*pi^22)/40343570821929375 + x2*((174611*pi^20)/183759535834875000 + x2*((-43867*pi^18)/196455649219830000 + x2*((3617*pi^16)/118168811560800000 + x2*(-1/364105581840000*pi^14 + x2*((691*pi^12)/3976032953692800000 + x2*(-1/122339475498240000*pi^10 + x2*(pi^8/3361246195507200000 + x2*(-1/114954619886346240000*pi^6 + x2*(pi^4/4598184795453849600000 + x2*(-1/155112100433309859840000*pi^2 + (pi/1240896803466478878720000 - 1/31022420086661971968000000*x)*x)))))))))))))

    sgn*sl_series(n, x)
end

function sl(n::Integer, z::Complex)
    n < 1  && throw(DomainError(n, "sl(n,z) undefined for n < 1"))

    if n == 1 && z == zero(z)
        zero(z)
    else
        eiz = exp(im*z)

        if iseven(n)
            (PolyLog.li(n, inv(eiz)) + PolyLog.li(n, eiz))/2
        else
            (PolyLog.li(n, inv(eiz)) - PolyLog.li(n, eiz))*im/2
        end
    end
end
