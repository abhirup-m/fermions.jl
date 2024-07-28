#### Iterative diagonalisation solution of the 1D Hubbard model ####
using Plots, Measures, Dates, ProgressMeter
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
        ("nn", [upPosition, upPosition + 1], U),
    ]
    return [hopping; correlation]
end

const t = 1.0
const retainSizeY = 50
const retainSizeX = 1e-2

function main(numSteps, U)
    initBasis = BasisStates(4; totOccCriteria=(x, N) -> div(N, 2) - 1 ≤ sum(x) ≤ div(N, 2) + 1, magzCriteria=x -> -2 ≤ sum(x[1:2:end] - x[2:2:end]) ≤ 2)
    hamiltonianFamily = [dimerHamiltonian(U, t)]
    numStatesFamily = Int64[2]
    for i in 1:numSteps
        push!(hamiltonianFamily, dimerAdditionalHamiltonian(U, t, 2 + i)
        )
        push!(numStatesFamily, 2 + i)
    end

    eigValData, eigVecData = IterDiag(
        hamiltonianFamily,
        initBasis,
        numStatesFamily,
        retainSizeY,
        retainSizeX;
        showProgress=true,
        totOccCriteria=(o, N) -> N - 1 ≤ o ≤ N + 1,
        magzCriteria=m -> -2 ≤ m ≤ 2
    )
    halfFillingSpectrum = [(E, X) for (E, X, numSites) in zip(eigValData, eigVecData, numStatesFamily)
                           if (numSites, 0) ∈ keys(E)]
    groundStateInfoAll = [(E[(numSites, 0)][1], X[(numSites, 0)][1]) for (E, X, numSites) in zip(eigValData, eigVecData, numStatesFamily)
                           if (numSites, 0) ∈ keys(E)]
    specfuncFlow, freqArray = IterSpecFunc(groundStateInfoAll, halfFillingSpectrum, (("-", [1], 1.0), ("+", [1], 1.0)), 5.5, 0.01)
    plots = []
    for y in specfuncFlow
        push!(plots, plot(freqArray, y))
    end
    p = plot(plots..., thickness_scaling=1.4, linewidth=2, layout=(length(plots), 1), leftmargin=4mm, size=(400, 300 * length(plots)), legend=false)
    display(p)
    #=savefig(p, "specfunc.pdf")=#
end

main(6, 2.0)
