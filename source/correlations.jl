using LinearAlgebra
include("fermionise.jl")

function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    minimumEnergies = minimum.(values(eigvals))
    minimumBlock = collect(keys(eigvals))[argmin(minimumEnergies)]
    minimumIndex = argmin(eigvals[minimumBlock])
    gstate = eigvecs[minimumBlock][minimumIndex]
    operatorMatrix = generalOperatorMatrix(basisStates, operatorList)
    return gstate' * operatorMatrix[minimumBlock] * gstate
end
