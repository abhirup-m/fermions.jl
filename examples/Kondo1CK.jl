using fermions, Plots, Measures
include("../src/iterDiag.jl")

realSites = 2
initSites = 1
maxSize = 400
V = 1.1
t = 1.
U = 0.

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

savePaths = IterDiag(hamFlow, maxSize);
f = deserialize(savePaths[end])
basicMats = f["operators"]
rotation = f["basis"]
eigVals = f["eigVals"]
display(eigVals[3][1:5])

basis = BasisStates(2 * (1 + realSites); totOccReq=[3])
fullHam = [
           [("n", [1], -U/2)];
           [("n", [2], -U/2)];
           [("nn", [1, 2], U)];
           [("+-",  [i, i+2], -V) for i in 1:2*initSites];
           [("+-",  [i+2, i], -V) for i in 1:2*initSites];
          ]
for site in initSites+1:realSites
    newTerm = [
               ("+-",  [2 * site - 1, 2 * site + 1], -t),
               ("+-",  [2 * site + 1, 2 * site - 1], -t),
               ("+-",  [2 * site, 2 * site + 2], -t),
               ("+-",  [2 * site + 2, 2 * site], -t),
              ]
    append!(fullHam, newTerm)
end
fullMatrix = OperatorMatrix(basis, fullHam)
E, X = eigen(fullMatrix)
display(E[1:5])

# display([IterCorrelation(path, [("+-+-",  [1, 2, 4, 3], 2.)], 1 + i) for (i, path) in enumerate(savePaths)])
