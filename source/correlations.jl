using LinearAlgebra
include("fermionise.jl")

function gstateCorrelation(basisStates, eigvals, eigvecs, operatorList)
    minimumBlock = collect(keys(eigvals))[argmin(minimum.(values(eigvals)))]
    minimumIndex = argmin(eigvals[minimumBlock])
    gstate = eigvecs[minimumBlock][minimumIndex]
    operatorMatrix = generalOperatorMatrix(basisStates, operatorList)
    return gstate' * operatorMatrix[minimumBlock] * gstate
end
