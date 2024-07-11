#### Iterative diagonalisation solution of the 1D Hubbard model ####
using Plots, Measures, BenchmarkTools
include("../src/base.jl")
include("../src/correlations.jl")
include("../src/eigen.jl")
include("../src/iterativeDiag.jl")

function dimerHamiltonian(U, t)
    hopping = [
        ("+-", [1, 3], -t),
        ("+-", [2, 4], -t),
        ("+-", [3, 1], -t),
        ("+-", [4, 2], -t),
    ]
    correlation = [
        ("n", [1], -U / 2),
        ("n", [2], -U / 2),
        ("nn", [1, 2], U),
        ("n", [3], -U / 2),
        ("n", [4], -U / 2),
        ("nn", [3, 4], U),
    ]
    return [hopping; correlation]
end


function dimerAdditionalHamiltonian(U, t, index)
    upPosition = 2 * index - 1
    hopping = [
        ("+-", [upPosition, upPosition - 2], -t),
        ("+-", [upPosition + 1, upPosition - 1], -t),
        ("+-", [upPosition - 2, upPosition], -t),
        ("+-", [upPosition - 1, upPosition + 1], -t),
    ]

    correlation = [
        ("n", [upPosition], -U / 2),
        ("n", [upPosition + 1], -U / 2),
        ("nn", [upPosition, upPosition + 1], -U),
    ]
    return [hopping; correlation]
end

function main(numSteps, U)
    t = 1.0
    initBasis = BasisStates(4)#; occCriteria=x -> x ∈ [1, 2, 3], magzCriteria = x -> x ∈ [-1, 0, 1])#; magzCriteria=magzCriteria, occCriteria=occCriteria)
    retainSize = 50
    hamiltonianFamily = [dimerHamiltonian(U, t)]
    numStatesFamily = Int64[2]
    for i in 1:numSteps
        push!(hamiltonianFamily, dimerAdditionalHamiltonian(U, t, 2 + i)
        )
        push!(numStatesFamily, 2 + i)
    end

    @time IterDiag(
        hamiltonianFamily,
        initBasis,
        numStatesFamily,
        retainSize;
        #occCriteria= (x,y) -> x ∈ [y-1, y, y+1], magzCriteria = x -> x ∈ [-1, 0, 1]
    )
    return
end

@btime main(7, 0.0);
