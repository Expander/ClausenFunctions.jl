module ClausenFunctionsForwardDiffExt

import ClausenFunctions
import ForwardDiff

# d/dx Cl₂(x) = Cl₁(x) = -log|2 sin(x/2)|
function ClausenFunctions.cl2(d::ForwardDiff.Dual{T}) where T
    x = ForwardDiff.value(d)
    ForwardDiff.Dual{T}(ClausenFunctions.cl2(x), ClausenFunctions.cl1(x)*ForwardDiff.partials(d))
end

end
