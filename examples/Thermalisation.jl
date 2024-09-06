using Plots, ProgressMeter, LinearAlgebra, ProgressMeter, fermions, Measures
theme(:dark)

include("../src/thermalisation.jl")

onsite = 2.

function HamiltonianAndState(numSites, onsite)
    @assert numSites % 2 == 0
    hamiltonian = [
                   [("n", [i], onsite) for i in 1:numSites];
                   [("+-", [i, i+1], -1.0) for i in 1:(numSites-1)];
                   [("+-", [i+1, i], -1.0) for i in 1:(numSites-1)]
                  ]
    initState = Dict(Bool.(vcat(repeat([1], div(numSites, 2)), 
                                repeat([0], div(numSites, 2)))
                          ) => 1.0
                    )
    return hamiltonian, initState
end

localOperator = [("n", [1], 1.0)]
localOperator2 = [("n", [2], 1.0)]
hopOperator = [("+-", [1, 2], 1.0), ("+-", [2, 1], 1.0)]
p1 = [] 
p2 = []
for numSites in [2, 6, 10]
    hamiltonian, initState = HamiltonianAndState(numSites, onsite)
    basisStates = BasisStates(numSites)
    @time occTimeEvol, timeVals, operatorMatrixTimeVol, hamiltonianMatrix, basisStates = OperatorTimeEvol(localOperator, hamiltonian, initState, basisStates, 0.1, 500; occupancy=div(numSites, 2))
    staticOperatorMatrix = OperatorMatrix(basisStates, localOperator2)
    @time otoc = OTOC(operatorMatrixTimeVol, staticOperatorMatrix, hamiltonianMatrix)
    push!(p1, plot(timeVals, occTimeEvol, label="\$N_s=$(numSites)\$"))
    push!(p2, plot(timeVals, otoc, label="\$N_s=$(numSites)\$"))
end
p1 = plot(p1..., layout=(length(p1), 1), linewidth=2, thickness_scaling=2.1, size=(800, 1500), leftmargin=-6mm, bottommargin=-4mm, ylabel="\$\\langle O(t)\\rangle\$", xlabel="time \$(t)\$")
p2 = plot(p2..., layout=(length(p2), 1), linewidth=2, thickness_scaling=2.1, size=(800, 1500), leftmargin=-6mm, bottommargin=-4mm, ylabel="\$\\langle [n_1(t), n_2(0)]^2\\rangle\$", xlabel="time \$(t)\$")
savefig(p1, "occupancyTime.pdf")
savefig(p2, "otoc.pdf")

