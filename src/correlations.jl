function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, correlationOperatorList::Vector{Dict{Tuple{String,Vector{Int64}},Float64}})
    gsEnergy, groundStates, blocks = getGstate(eigvals, eigvecs)
    correlationResults = []
    for correlationOperator in correlationOperatorList
        correlationMatrix = generalOperatorMatrix(basisStates, correlationOperator)
        correlationArray = [state' * correlationMatrix[block] * state for (state, block) in zip(groundStates, blocks)]
        push!(correlationResults, sum(abs.(correlationArray ./ length(correlationArray))))
    end
    return correlationResults
end
