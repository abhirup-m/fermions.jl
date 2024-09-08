using Plots, ProgressMeter, LinearAlgebra, ProgressMeter, fermions, Measures
theme(:dark)

include("../src/thermalisation.jl")

t = 1.
U = 2.

function HamiltonianAndState(numSites)
    @assert (2 + numSites) % 2 == 0
    hamiltonian = [
                   [("nn", [1, 2], U)];
                   [("+-", [i, i+1], -t) for i in 1:numSites+1];
                   [("+-", [i+1, i], -t) for i in 1:numSites+1];
                  ]
    initState = Dict(Bool.([ifelse(i % 2 == 0, 1, 0) for i in 1:2+numSites]) => 1.0
                    )
    return hamiltonian, initState
end

impOccOperator = [("n", [2], 1.0)]#, ("n", [1], 1.0)]
bathLocalOccOperator = [("n", [3], 1.0)]#, ("n", [4], 1.0)]
hopOperator = [("+-", [1, 2], 1.0), ("+-", [2, 1], 1.0)]
p1 = [] 
p2 = []
for numSites in [2, 4, 6, 12]
    hamiltonian, initState = HamiltonianAndState(numSites)
    basisStates = BasisStates(2 + numSites, totOccReq=div(2 + numSites, 2))
    hamiltonianMatrix = OperatorMatrix(basisStates, hamiltonian)
    @time occTimeEvol, timeVals, operatorMatrixTimeVol = OperatorTimeEvol(impOccOperator, hamiltonian, initState, basisStates, 0.1, 100)
    staticOperatorMatrix = OperatorMatrix(basisStates, bathLocalOccOperator)
    @time otoc = OTOC(operatorMatrixTimeVol, staticOperatorMatrix, hamiltonianMatrix)
    push!(p1, plot(timeVals, occTimeEvol, label="\$N_s=$(numSites)\$"))
    push!(p2, plot(timeVals, otoc, label="\$N_s=$(numSites)\$"))
end
p1 = plot(p1..., layout=(length(p1), 1), linewidth=2, thickness_scaling=2.1, size=(800, 1500), leftmargin=-6mm, bottommargin=-4mm, ylabel="\$\\langle O(t)\\rangle\$", xlabel="time \$(t)\$")
p2 = plot(p2..., layout=(length(p2), 1), linewidth=2, thickness_scaling=2.1, size=(800, 1500), leftmargin=-6mm, bottommargin=-4mm, ylabel="\$\\langle [n_1(t), n_2(0)]^2\\rangle\$", xlabel="time \$(t)\$")
savefig(p1, "occupancyTime.pdf")
savefig(p2, "otoc.pdf")
