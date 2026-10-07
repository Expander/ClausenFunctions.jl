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

end
