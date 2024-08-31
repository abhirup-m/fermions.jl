using fermions, Plots, Measures
include("../src/iterDiag.jl")

realSites = 5
initSites = 1
maxSize = 500
V = 1.1
t = 1.
U = 1.

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
display(f["eigVals"][1+realSites][1:5])
corrdef = [("+-+-", [1,2,4,3], 1.0)]
println(IterCorrelation(savePaths[end], corrdef, 1 + realSites))

#=basis = BasisStates(2 * (1 + realSites); totOccReq=[1+realSites])=#
#=fullHam = [=#
#=           [("n", [1], -U/2)];=#
#=           [("n", [2], -U/2)];=#
#=           [("nn", [1, 2], U)];=#
#=           [("+-",  [i, i+2], -V) for i in 1:2*initSites];=#
#=           [("+-",  [i+2, i], -V) for i in 1:2*initSites];=#
#=          ]=#
#=for site in initSites+1:realSites=#
#=    newTerm = [=#
#=               ("+-",  [2 * site - 1, 2 * site + 1], -t),=#
#=               ("+-",  [2 * site + 1, 2 * site - 1], -t),=#
#=               ("+-",  [2 * site, 2 * site + 2], -t),=#
#=               ("+-",  [2 * site + 2, 2 * site], -t),=#
#=              ]=#
#=    append!(fullHam, newTerm)=#
#=end=#
#=fullMatrix = OperatorMatrix(basis, fullHam)=#
#=E, X = eigen(fullMatrix)=#
#=display(E[1:5])=#
#=display(X[:, 1]' * OperatorMatrix(basis, corrdef) * X[:, 1])=#
