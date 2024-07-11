include("../src/base.jl")
include("../src/correlations.jl")
include("../src/iterativeDiag.jl")
include("testing_helpers.jl")

using Test
using LinearAlgebra

include("base_tests.jl")
# include("eigen_tests.jl")
# include("correlation_tests.jl")
# include("model_tests.jl")
# include("iterDiag_tests.jl")
