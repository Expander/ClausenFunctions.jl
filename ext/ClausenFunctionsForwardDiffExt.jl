module ClausenFunctionsForwardDiffExt

import ClausenFunctions
import ForwardDiff

# d/dx Cl₁(x) = -cot(x/2)/2
function ClausenFunctions.cl1(d::ForwardDiff.Dual{T}) where T
    cl0(x) = cot(x/2)/2
    x = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl1(x), -cl0(x)*ForwardDiff.partials(d))
end

# d/dx Cl₂(x) = Cl₁(x) = -log|2 sin(x/2)|
function ClausenFunctions.cl2(d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl2(x), ClausenFunctions.cl1(x)*ForwardDiff.partials(d))
end

# d/dx Cl3(x) = -Cl2(x)
function ClausenFunctions.cl3(d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl3(x), -ClausenFunctions.cl2(x)*ForwardDiff.partials(d))
end

# d/dx Cl4(x) = Cl3(x)
function ClausenFunctions.cl4(d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl4(x), ClausenFunctions.cl3(x)*ForwardDiff.partials(d))
end

# d/dx Cl5(x) = -Cl4(x)
function ClausenFunctions.cl5(d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl5(x), -ClausenFunctions.cl4(x)*ForwardDiff.partials(d))
end

# d/dx Cl6(x) = Cl5(x)
function ClausenFunctions.cl6(d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl6(x), ClausenFunctions.cl5(x)*ForwardDiff.partials(d))
end

# d/dx Cl_{2n+2}(x) = Cl_{2n+1}(x), d/dx Cl_{2n+1}(x) = -Cl_{2n}(x)
function ClausenFunctions.cl(n::Integer, d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    sgn = iseven(n) ? 1 : -1
    ForwardDiff.Dual{T}(ClausenFunctions.cl(n,x), flipsign(ClausenFunctions.cl(n-1,x),sgn)*ForwardDiff.partials(d))
end

# d/dx Sl_{2n+2}(x) = -Sl_{2n+1}(x), d/dx Sl_{2n+1}(x) = Sl_{2n}(x)
function ClausenFunctions.sl(n::Integer, d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    sgn = iseven(n) ? -1 : 1
    ForwardDiff.Dual{T}(ClausenFunctions.sl(n,x), flipsign(ClausenFunctions.sl(n-1,x),sgn)*ForwardDiff.partials(d))
end

end
