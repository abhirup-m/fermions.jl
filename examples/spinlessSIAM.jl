using fermions, Plots, Measures, BenchmarkTools
theme(:dark)
include("../src/iterDiag.jl")

realSites = 18
initSites = 3
V = 1.0
t = 1.0
occupancy = 2

function getHamFlow(initSites, realSites)
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[]
    initHam = [
               [("+-",  [i, i+1], -V) for i in 1:initSites];
               [("+-",  [i+1, i], -V) for i in 1:initSites];
              ]
    push!(hamFlow, initHam)
    for site in initSites+1:realSites
        newTerm = [
                   ("+-",  [site, site + 1], -t),
                   ("+-",  [site + 1, site], -t),
                  ]
        push!(hamFlow, newTerm)
    end
    return hamFlow
end

hamFlow = getHamFlow(initSites, realSites)
p1 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="spin-flip corr.")
p2 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="energy per site")
corrdef = [("+-", [1,2], 1.0), ("+-", [2,1], 1.0)]
for maxSize in [30, 100, 800]
    savePaths = IterDiag(hamFlow, maxSize; symmetries=Char['N'], occReq=(x,N)->ifelse(0 ≤ x ≤ 4, true, false))
    corrFlow, energypersite = IterCorrelation(savePaths, corrdef; occupancy=occupancy)
    scatter!(p1, initSites:realSites, corrFlow, label="M=$(maxSize)")
    scatter!(p2, initSites:realSites, energypersite, label="M=$(maxSize)")
end

corrExact = []
energyExact = []
@showprogress for (i, num) in enumerate(initSites:realSites)
    basis = BasisStates(1 + num; totOccReq=[occupancy])
    fullHam = vcat(hamFlow[1:i]...)
    fullMatrix = OperatorMatrix(basis, fullHam)
    F = eigen(fullMatrix)
    push!(corrExact, F.vectors[:, 1]' * OperatorMatrix(basis, corrdef) * F.vectors[:, 1])
    push!(energyExact, minimum(F.values)/(1 + num))
end

plot!(p1, initSites:realSites, corrExact, label="exact", linewidth=2)
plot!(p2, initSites:realSites, energyExact, label="exact", linewidth=2)
savefig(p1, "spinflipcomparison.pdf")
savefig(p2, "energy.pdf")
