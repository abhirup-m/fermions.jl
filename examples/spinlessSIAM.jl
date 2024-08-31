using fermions, Plots, Measures
include("../src/iterDiag.jl")

maxSize = 100
V = 1.
t = 2.
E = 1.
occupancy = 2

hamFlow = [
           [
            ("+-",  [1, 2], -V),
            ("+-",  [2, 1], -V),
           ],
           [
            ("+-",  [2, 3], -t),
            ("+-",  [3, 2], -t),
           ],
          ]
savePaths = IterDiag(hamFlow, maxSize);
f = deserialize(savePaths[end])
basicMats = f["operators"]
rotation = f["basis"]
eigVals = f["eigVals"]
display(eigVals)

basis = BasisStates(3)
fullHam = vcat(hamFlow...)
fullMatrix = OperatorMatrix(basis, fullHam)
E, X = eigen(fullMatrix)
display(E)
