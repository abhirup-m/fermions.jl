using Test, LinearAlgebra, ProgressMeter, Serialization
using fermions

include("testing_helpers.jl")
include("base_tests.jl")
include("eigen_tests.jl")
include("correlation_tests.jl")
include("eigenstateRG_tests.jl")
include("iterDiag_tests.jl")
#=include("modelHamiltonians.jl")=#
