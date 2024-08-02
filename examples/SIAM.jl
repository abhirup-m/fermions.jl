using fermions
include("../src/iterDiag.jl")

realSites = 20
maxSize = 400
U = 1.
t = 1.
hamFlow = [
           [
           ("n", [1], -U/2),
          ]
          ]
for site in 1:realSites
    newTerm = [
               ("+-",  [1, 1 + site], -t),
               ("+-",  [1 + site, 1], -t),
              ]
    push!(hamFlow, newTerm)
end

initBasis = BasisStates(1)
numSitesFlow = collect(1:1:1+realSites)
r = IterDiag(hamFlow, maxSize);
