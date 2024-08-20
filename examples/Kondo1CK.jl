using fermions, Plots, Measures
include("../src/iterDiag.jl")

realSites = 10
maxSize = 200
J = 1.
t = 0.1

hamFlow = [
           [
           ("nn",  [1, 3], J/4), # z-uu
           ("nn",  [1, 4], -J/4), # z-ud
           ("nn",  [2, 3], -J/4), # z-du
           ("nn",  [2, 4], J/4), # z-dd
           ("+-+-",  [1, 2, 4, 3], J/2), # +-
           ("+-+-",  [2, 1, 3, 4], J/2), # -+
          ]
          ]
for site in 2:realSites
    newTerm = [
               ("+-",  [2 * site - 1, 2 * site + 1], -t),
               ("+-",  [2 * site + 1, 2 * site - 1], -t),
               ("+-",  [2 * site, 2 * site + 2], -t),
               ("+-",  [2 * site + 2, 2 * site], -t),
              ]
    push!(hamFlow, newTerm)
end

initBasis = BasisStates(1)
numSitesFlow = collect(1:1:1+realSites)
savePaths = IterDiag(hamFlow, maxSize);
totOcc = Tuple{String, Vector{Int64}, Float64}[] #[("n", [i], 1.0) for i in 1:2*(1+realSites)]
totSz = Tuple{String, Vector{Int64}, Float64}[] #[("n", [i], (-1.)^i) for i in 1:2*(1+realSites)]
for i in 1:2*(1+realSites)
    push!(totOcc, ("n", [i], 1.0))
    push!(totSz, ("n", [i], (-1.)^i))
end
spinOperators = [[("+-+-",  [1, 2, 2*i+2, 2*i+1], 1.)] for i in 1:realSites]
println([IterCorrelation(path, totOcc[1:i]) for (i, path) in enumerate(savePaths)])
println([IterCorrelation(path, totSz[1:i]) for (i, path) in enumerate(savePaths)])
spinCorr = [IterCorrelation(savePaths[end], operator) for operator in spinOperators]
println(spinCorr)
