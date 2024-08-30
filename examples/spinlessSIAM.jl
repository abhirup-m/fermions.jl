using fermions, Plots, Measures
include("../src/iterDiag.jl")

realSites = 2
initSites = 1
maxSize = 100
V = 1.
t = V

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
               ("+-",  [site + 1, site + 2], -t),
               ("+-",  [site + 2, site + 1], -t),
              ]
    push!(hamFlow, newTerm)
end

savePaths = IterDiag(hamFlow, maxSize);
println("------------")
display([IterCorrelation(path, [("+-",  [1, 2], 1.)], i) for (i, path) in enumerate(savePaths)])
