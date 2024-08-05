using fermions, Plots, Measures
include("../src/iterDiag.jl")

realSites = 20
maxSize = 200
U = 8.
t = 1.1
# hamFlow = [[("n", [1], -U/2)]]
# for site in 1:realSites
#     newTerm = [
#                ("+-", [site, 1 + site], -t),
#                ("+-", [1 + site, site], -t),
#               ]
#     push!(hamFlow, newTerm)
# end

hamFlow = [
           [
           ("n", [1], -U/2),
           ("n", [2], -U/2),
           ("nn", [1, 2], U),
           ("+-",  [1, 3], -t),
           ("+-",  [3, 1], -t),
           ("+-",  [2, 4], -t),
           ("+-",  [4, 2], -t),
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

broadening = 10^-1.5
freqArray = collect(range(-10, stop=10, step=0.01))
probe = [("-", [1], 1.0), ("-", [2], 1.0)]
probeDag = [("+", [1], 1.0), ("+", [2], 1.0)]
specfuncFlow, freqArray = IterSpecFunc(savePaths[1:2:end], probe, probeDag, freqArray, broadening)

for (i, SF) in enumerate(specfuncFlow)
    fig = plot(freqArray, SF, size=(350, 200), thickness_scaling=1.3)
    savefig(fig, "specfunc_$(i).pdf")
end
run(`pdfunite $(["specfunc_$(i).pdf" for i in eachindex(specfuncFlow)]) specfunc.pdf`)
run(`rm $(["specfunc_$(i).pdf" for i in eachindex(specfuncFlow)])`)
