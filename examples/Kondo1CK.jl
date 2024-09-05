using fermions, Plots, Measures
theme(:dark)
include("../src/iterDiag.jl")

realSites = 7
initSites = 1
maxSize = 50
V = 1.1
t = 1.1
U = 1.

function getHamFlow(initSites, realSites)
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[]
    initHam = [
               [("n", [1], -U/2)];
               [("n", [2], -U/2)];
               [("nn", [1, 2], U)];
               [("+-",  [i, i+2], -V) for i in 1:2*initSites];
               [("+-",  [i+2, i], -V) for i in 1:2*initSites];
              ]
    push!(hamFlow, initHam)
    for site in initSites+1:realSites
        newTerm = [
                   ("+-",  [2 * site - 1, 2 * site + 1], -t),
                   ("+-",  [2 * site + 1, 2 * site - 1], -t),
                   ("+-",  [2 * site, 2 * site + 2], -t),
                   ("+-",  [2 * site + 2, 2 * site], -t),
                  ]
        push!(hamFlow, newTerm)
    end
    return hamFlow
end

hamFlow = getHamFlow(initSites, realSites)
savePaths = IterDiag(hamFlow, maxSize);
corrdef = [("+-+-", [1,2,4,3], 1.0)]
corrFlow = IterCorrelation(savePaths, corrdef)

corrExact = []
for num in initSites:realSites
    basis = BasisStates(2 * (1 + num); totOccReq=[1 + num], magzReq = [-1, 0, 1])
    fullHam = vcat(hamFlow[initSites:num]...)
    fullMatrix = OperatorMatrix(basis, fullHam)
    F = eigen(fullMatrix)
    push!(corrExact, F.vectors[:, 1]' * OperatorMatrix(basis, corrdef) * F.vectors[:, 1])
end

p = scatter(initSites:realSites, corrFlow, thickness_scaling=2, leftmargin=-12mm, bottommargin=-6mm, label="approx.", xlabel="sites", ylabel="spin-flip corr.")
plot!(initSites:realSites, corrExact, thickness_scaling=2, label="exact")
savefig(p, "corr.pdf")

# realSites = 30
# hamFlow = getHamFlow(initSites, realSites)
# savePaths = IterDiag(hamFlow, 200);
# energypersite = []
# for (step, savePath) in enumerate(savePaths[1:end-1])
#     f = deserialize(savePath)
#     basis = f["basis"]
#     eigVals = f["eigVals"]
#     push!(energypersite, eigVals[1] / (step + 1))
# end
# p = scatter(initSites:realSites, energypersite, thickness_scaling=1.6, leftmargin=-12mm, bottommargin=-6mm, xlabel="sites", ylabel="energy per site")
# savefig(p, "energy.pdf")
# 
# corrdef = [("+-+-", [1,2,4,3], 1.0)]
# corrFlow = IterCorrelation(savePaths, corrdef)
# p = scatter(initSites:realSites, corrFlow, thickness_scaling=1.6, leftmargin=-12mm, bottommargin=-6mm, xlabel="sites", ylabel="spin-flip corr.")
# savefig(p, "corrApprox.pdf")
