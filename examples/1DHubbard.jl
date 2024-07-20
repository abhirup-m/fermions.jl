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
        ("nn", [upPosition, upPosition + 1], U),
    ]
    return [hopping; correlation]
end

function main(numSteps, U)
    t = 1.0
    initBasis = BasisStates(4; totOccCriteria=(x, N) -> div(N, 2) - 1 ≤ sum(x) ≤ div(N, 2) + 1, magzCriteria=x -> -2 ≤ sum(x[1:2:end] - x[2:2:end]) ≤ 2)
    retainSize = 200
    hamiltonianFamily = [dimerHamiltonian(U, t)]
    numStatesFamily = Int64[2]
    for i in 1:numSteps
        push!(hamiltonianFamily, dimerAdditionalHamiltonian(U, t, 2 + i)
        )
        push!(numStatesFamily, 2 + i)
    end

    spectrumFlow = IterDiag(
        hamiltonianFamily,
        initBasis,
        numStatesFamily,
        retainSize;
        totOccCriteria=(o, N) -> N - 1 ≤ o ≤ N + 1,
        magzCriteria=m -> -2 ≤ m ≤ 2
    )

    freqArray = collect(range(-10.0, stop=10.0, step=0.01))
    specfunc = 0 .* freqArray
    plots = []
    for (numSites, spectrum) in zip(numStatesFamily, spectrumFlow)
        probe, probeDag = (("-", [numSites], 1.0), ("+", [numSites], 1.0))
        if (numSites, 0) ∉ keys(spectrum)
            continue
        end
        eigVals = [spectrum[(numSites, 0)][1]; vcat([v[1] for (k,v) in spectrum if k ≠ (numSites, 0)]...)]
        eigVecs = [spectrum[(numSites, 0)][2]; vcat([v[2] for (k,v) in spectrum if k ≠ (numSites, 0)]...)]
        specfunc .+= SpecFunc(eigVals, eigVecs, probe, probeDag, freqArray, 1e-2)
        push!(plots, plot(freqArray, specfunc))
    end
    specfunc ./= sum(specfunc) * abs(freqArray[2] - freqArray[1])
    return plots
end

@time plots = main(10, 2.0);
p = plot(plots..., thickness_scaling=1.4, linewidth=2, layout=(length(plots), 1), leftmargin=4mm, size=(400, 300 * length(plots)), legend=false)
savefig(p, "specfunc.pdf")
