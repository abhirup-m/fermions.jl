using LinearAlgebra
include("fermionise.jl")

function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, correlationOperator::Dict{Tuple{String,Vector{Int64}},Float64})
    gsEnergy = minimum(round.(minimum.(values(eigvals)), digits=trunc(Int, -log10(TOLERANCE))))
    correlationMatrix = generalOperatorMatrix(basisStates, correlationOperator)
    correlations = [eigvecs[block][i]' * correlationMatrix[block] * eigvecs[block][i]
                    for (block, eigvalsBlock) in eigvals
                    for (i, _) in enumerate(eigvalsBlock[abs.(eigvalsBlock .- gsEnergy).<TOLERANCE])]
    return sum(abs.(correlations ./ length(correlations)))
end


function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, correlationOperatorList::Vector{Dict{Tuple{String,Vector{Int64}},Float64}})
    gsEnergy = minimum(round.(minimum.(values(eigvals)), digits=trunc(Int, -log10(TOLERANCE))))
    correlationResults = []
    for correlationOperator in correlationOperatorList
        correlationMatrix = generalOperatorMatrix(basisStates, correlationOperator)
        correlations = [eigvecs[block][i]' * correlationMatrix[block] * eigvecs[block][i]
                        for (block, eigvalsBlock) in eigvals
                        for (i, _) in enumerate(eigvalsBlock[abs.(eigvalsBlock .- gsEnergy).<TOLERANCE])]
        push!(correlationResults, sum(abs.(correlations ./ length(correlations))))
    end
    return correlationResults
end
