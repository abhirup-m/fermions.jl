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
    initBasis = BasisStates(4, [1, 2, 3], [-1, 0, 1], x -> true)
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
    groundStateInfoAll = [(numSites, 0) for (E, X, numSites) in zip(eigValData, eigVecData, numStatesFamily)
                           if (numSites, 0) ∈ keys(E)]

    freqArray = collect(range(-6, stop=6, step=0.01))
    probe = [("-", [1], 1.0)]
    probeDag = [("+", [1], 1.0)]
    specfuncFlow, freqArray = IterSpecFunc(groundStateInfoAll, halfFillingSpectrum, probe, probeDag, freqArray, 0.01, ['N', 'Z'])
    plots = []
    for y in specfuncFlow
        push!(plots, plot(freqArray, y))
    end
    p = plot(plots..., thickness_scaling=1.4, linewidth=2, layout=(length(plots), 1), leftmargin=4mm, size=(400, 300 * length(plots)), legend=false)
    display(p)
    #=savefig(p, "specfunc.pdf")=#
end

main(6, 2.0)
