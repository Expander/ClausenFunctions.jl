import SpecialFunctions

# (-1)^k B_{2k}/(2k)! where B_{2k} are the even Bernoulli numbers
# B(k) = 2*(-1)^(2*k + 1)*SpecialFunctions.zeta(2*k)/(2*pi)^(2*k)
const B = (-0.083333333333333333,-0.0013888888888888889,-0.000033068783068783069,-8.2671957671957672e-7,
           -2.0876756987868099e-8,-5.2841901386874932e-10,-1.3382536530684679e-11,-3.3896802963225829e-13,
           -8.5860620562778446e-15,-2.1748686985580619e-16,-5.5090028283602295e-18,-1.3954464685812523e-19,
           -3.5347070396294675e-21,-8.9535174270375469e-23,-2.2679524523376831e-24,-5.7447906688722024e-26,
           -1.4551724756148649e-27,-3.6859949406653102e-29,-9.3367342570950447e-31,-2.3650224157006299e-32,
           -5.9906717624821343e-34,-1.5174548844682903e-35,-3.8437581254541882e-37,-9.736353072646691e-39,
           -2.466247044200681e-40,-6.2470767418207437e-42,-1.5824030244644914e-43,-4.008273685948936e-45,
           -1.0153075855569556e-46,-2.5718041582418717e-48,-6.5144560352338149e-50,-1.6501309906896525e-51,
           -4.1798306285394759e-53,-1.0587634667702909e-54,-2.6818791912607707e-56,-6.7932793511074212e-58,
           -1.7207577616681405e-59,-4.3587303293488938e-61,-1.1040792903684667e-62,-2.7966655133781345e-64,
           -7.0840365016794702e-66,-1.7944074082892241e-67,-4.5452870636110961e-69,-1.1513346631982052e-70,
           -2.9163647710923614e-72,-7.3872382634973376e-74,-1.8712093117637953e-75,-4.7398285577617994e-77,
           -1.2006125993354507e-78,-3.0411872415142924e-80,-7.7034172747051063e-82,-1.9512983909098831e-83,
           -4.9426965651594615e-85,-1.2519996659171848e-86,-3.1713522017635155e-88,-8.0331289707353345e-90,
           -2.0348153391661466e-91,-5.1542474664474739e-93,-1.3055861352149467e-94,-3.3070883141750912e-96,
           -8.3769525600490913e-98,-2.1219068717497138e-99,-5.3748528956122803e-101,-1.3614661432172069e-102,
           -3.448634027993399e-104,-8.7354920416383551e-106,-2.2127259833925497e-107,-5.6049003928372242e-109,
           -1.4197378549991788e-110,-3.5962379982587627e-112,-9.1093772660782319e-114,-2.3074322171091233e-115,
           -5.844794085299002e-117,-1.4805036371705745e-118,-3.7501595226227197e-120,-9.4992650419929583e-122,
           -2.4061919444675199e-123,-6.0949553971026848e-125,-1.5438702377042471e-126,-3.9106689968592923e-128,
           -9.9058402898794298e-130,-2.5091786578563553e-131,-6.3558237896024598e-133,-1.609948973461625e-134,
           -4.0780483898724622e-136,-1.032981724531519e-137,-2.6165732752609202e-139,-6.6278575334086228e-141,
           -1.6788559257443673e-142,-4.2525917390343006e-144,-1.0771940713664573e-145,-2.728564458083958e-147,
           -6.9115345134370137e-149,-1.7507121442157682e-150,-4.4346056667239196e-152,-1.123298737848715e-153,
           -2.8453489425693971e-155,-7.2073530684151368e-157,-1.8256438595501418e-158,-4.6244099189746555e-160,
           -1.1713767165946985e-161,-2.9671318854112523e-163,-7.5158328663197353e-165,-1.9037827051837494e-166,
           -4.8223379271757316e-168,-1.2215124667619569e-169,-3.0941272241548313e-171,-7.8375158172837117e-173,
           -1.9852659485568249e-174,-5.0287373938151486e-176,-1.2737940624195893e-177,-3.2265580530233667e-179,
           -8.1729670255761088e-181,-2.0702367322529201e-182,-5.2439709032927841e-184,-1.3283133472690102e-185,
           -3.3646570148302941e-187,-8.5227757823275058e-189,-2.1588443254591854e-190,-5.4684165588767235e-192,
           -1.3851660959868732e-193,-3.5086667096656512e-195,-8.8875566007447618e-197,-2.2512443861893281e-198,
           -5.7024686469217722e-200,-1.4444521824735857e-201,-3.6588401210745441e-203,-9.2679502956336812e-205,
           -2.3475992347298971e-206,-5.9465383295169881e-208,-1.5062757553029781e-209,-3.8154410604763518e-211,
           -9.6646251091260105e-213,-2.4480781387902631e-214,-6.2010543667790175e-216,-1.5707454206813435e-217,
           -3.9787446306053875e-219,-1.0078277884588346e-220,-2.552857610857218e-222,-6.4664638700600961e-224,
           -1.637974433237253e-225,-4.1490377087871461e-227,-1.0509635290775169e-228,-2.6621217182765621e-230,
           -6.7432330873938832e-232,-1.7080808949773099e-233,-4.3266194508991167e-235,-1.09594550983765e-236,
           -2.7760624066064009e-238,-7.0318482225589324e-240,-1.7811879627573501e-241,-4.511801816901474e-243,
           -1.1428527511202687e-244,-2.8948798368101921e-246,-7.3328162891986569e-248,-1.8574240646335573e-249,
           -4.704910118860853e-251,-1.1917676554344843e-252,-3.0187827368818927e-254,-7.6466660014982318e-256,
           -1.9369231254735576e-257,-4.9062835924299293e-259,-1.2427761521749539e-260,-3.1479887685209103e-262,
           -7.9739487029830971e-264,-2.0198248022238287e-265,-5.1162759927867278e-267,-1.2959678485750701e-268,
           -3.2827249095010017e-270,-8.3152393350706909e-272,-2.1062747292467202e-273,-5.3352562160805552e-275,
           -1.3514361871210552e-276,-3.4232278524048296e-278,-8.6711374470768819e-280,-2.1964247741580717e-281,
           -5.5636089474762561e-283,-1.4092786097034883e-284,-3.5697444204246398e-286,-9.0422682494513886e-288,
           -2.2904333046148606e-289,-5.8017353369352214e-291,-1.4695967287946349e-292,-3.7225320009595016e-294,
           -9.4292837120924187e-296,-2.3884654665215509e-297,-6.0500537039203004e-299,-1.5324965059522865e-300,
           -3.8818589977708149e-302,-9.8328637096699497e-304,-2.490693474143869e-305,-6.3090002722625804e-307,
           -1.5980884379636898e-308,-4.0480053024903931e-310,-1.0253717215969654e-311,-2.5972969126396544e-313,
           -6.5790299364809836e-315,-1.6664877509565689e-316,-4.2212627863094242e-318,-1.069258354935559e-319,
           -2.7084630535382438e-321,-6.8606170609008831e-323)

# returns (x, sign) with x in [0,pi]
function range_reduce(n::Int64, x::Float64)
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

# returns N_n(x) from Eq.(2.11)
function ncal(n::Int64, x::Float64)
    sum = zero(x)
    xn = x^(n + 1)

    for k in one(n):length(B)
        xn *= x*x
        term = B[k]*xn/(2*k + n + 1)
        old_sum = sum
        sum += term
        sum == old_sum && break
    end

    (x^(n + 1)/(n + 1) + sum)/(n + 1)
end

# returns Cl(n, 0)
function cln0(n::Int64)::Float64
    if iseven(n)
        zero(Float64)
    else
        SpecialFunctions.zeta(n)
    end
end

# returns P_k(x)
function pcal(k::Int64, x::Float64)
    sum = zero(x)
    x2 = x*x

    for i in 3:2:k
        sum = x2*sum + (-1)^(fld(k - 1, 2.0) + fld(i - 1, 2.0))/factorial(k - i)*cln0(i)
    end

    if iseven(k)
        sum *= x
    end

    sum
end

# returns sum in Eq.(2.13)
function nsum(n::Int64, x::Float64)
    sum = zero(x)
    xn = one(x)

    for i in zero(n):(n - 2)
        sum += binomial(n - 2, i)*xn*ncal(n - 2 - i, x)
        xn *= -x
    end

    sum
end

"""
    cl(n::Int64, x::Float64)::Float64

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
function cl(n::Int64, x::Float64)::Float64
    if n < 2
        throw("cl(n,x) undefined for n < 2")
    end

    (x, sgn) = range_reduce(n, x)

    if x == zero(x)
        return cln0(n)
    end

    if iseven(n) && x == pi
        return zero(x)
    end

    fn2 = factorial(n - 2)

    # Eq.(2.13)
    sgn*((-1)^fld(n + 1, 2.0)*x^(n - 1)/(fn2*(n - 1))*log(2*sin(x/2))
         + (-1)^(fld(n, 2.0) + 1)/fn2*nsum(n, x) + pcal(n, x))
end
