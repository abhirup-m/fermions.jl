using fermions, Plots
include("../src/iterDiag.jl")

realSites = 10
maxSize = 200
U = 0.
t = 1.
hamFlow = [[("n", [1], -U/2)]]
for site in 1:realSites
    newTerm = [
               ("+-", [site, 1 + site], -t),
               ("+-", [1 + site, site], -t),
              ]
    push!(hamFlow, newTerm)
end

# hamFlow = [
#            [
#            ("n", [1], -U/2),
#            ("n", [2], -U/2),
#            ("nn", [1, 2], U),
#           ]
#           ]
# for site in 1:realSites
#     newTerm = [
#                ("+-",  [2 * site - 1, 2 * site + 1], -t),
#                ("+-",  [2 * site + 1, 2 * site - 1], -t),
#                ("+-",  [2 * site, 2 * site + 2], -t),
#                ("+-",  [2 * site + 2, 2 * site], -t),
#               ]
#     push!(hamFlow, newTerm)
# end

initBasis = BasisStates(1)
numSitesFlow = collect(1:1:1+realSites)
savePaths = IterDiag(hamFlow, maxSize);

broadening = 1e-3
freqArray = collect(range(-2.5, stop=2.5, step=0.01))
probe = [("-", [1], 1.0)]
probeDag = [("+", [1], 1.0)]
specfuncFlow, freqArray = IterSpecFunc(savePaths, probe, probeDag, freqArray, broadening)

plots = [plot(freqArray, specfuncFlow[i]) for i in eachindex(specfuncFlow)]
p = plot(plots..., layout=(length(plots), 1), size=(300, 200 * length(plots)))
savefig(p, "specfunc.pdf")
