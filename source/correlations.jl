using LinearAlgebra
include("fermionise.jl")

function gstateCorrelation(basisStates, eigvals, eigvecs, operatorList)
    minimumEnergies = minimum.(values(eigvals))
    minimumBlock = collect(keys(eigvals))[argmin(minimumEnergies)]
    minimumIndex = argmin(eigvals[minimumBlock])
    gstate = eigvecs[minimumBlock][minimumIndex]
    operatorMatrix = generalOperatorMatrix(basisStates, operatorList)
    return gstate' * operatorMatrix[minimumBlock] * gstate
end
