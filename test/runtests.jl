using Test
using LinearAlgebra, ProgressMeter

include("../src/base.jl")
include("../src/eigen.jl")
include("../src/correlations.jl")
include("../src/eigenstateRG.jl")
include("../src/reverseUnitaries.jl")
include("../src/stateExpansion.jl")
include("../src/iterativeDiag.jl")
include("testing_helpers.jl")

include("base_tests.jl")
include("eigen_tests.jl")
include("correlation_tests.jl")
include("eigenstateRG_tests.jl")
# include("iterDiag_tests.jl")
