using fermions, Plots, Measures
include("../src/iterDiag.jl")

realSites = 3
initSites = 1
maxSize = 100
V = 1.
t = V / 10
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
println("------------")
display([IterCorrelation(path, [("+-+-",  [1, 2, 4, 3], 2.)], 1 + i) for (i, path) in enumerate(savePaths)])
