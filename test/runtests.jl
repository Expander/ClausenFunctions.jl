using Test
using DelimitedFiles
import ForwardDiff
import ClausenFunctions

include("Cl2_reference.jl")

include("Cl1.jl")
include("Cl2.jl")
include("Cl3.jl")
include("Cl4.jl")
include("Cl5.jl")
include("Cl6.jl")
include("Cl.jl")
include("Missing.jl")
include("Sl.jl")
include("range_reduction.jl")
include("TypeStability.jl")
