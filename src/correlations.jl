function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, correlationOperator::Dict{Tuple{String,Vector{Int64}},Float64})
    gsEnergy, groundStates, blocks = getGstate(eigvals, eigvecs)
    correlationMatrix = generalOperatorMatrix(basisStates, correlationOperator)
    correlationArray = [state' * correlationMatrix[block] * state for (state, block) in zip(groundStates, blocks)]
    correlationResult = sum(abs.(correlationArray) ./ length(correlationArray))
    return correlationResult
end


function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, correlationOperatorList::Vector{Dict{Tuple{String,Vector{Int64}},Float64}})
    gsEnergy, groundStates, blocks = getGstate(eigvals, eigvecs)
    correlationResults = zeros(length(correlationOperatorList))
    Threads.@threads for (i, correlationOperator) in collect(enumerate(correlationOperatorList))
        correlationMatrix = generalOperatorMatrix(basisStates, correlationOperator)
        correlationArray = [state' * correlationMatrix[block] * state for (state, block) in zip(groundStates, blocks)]
        correlationResults[i] = sum(abs.(correlationArray) ./ length(correlationArray))
    end
    return correlationResults
end


function vnEntropy(basisStates::Vector{BitArray}, groundState::Vector{Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[])
    reducedDMatrix = reducedDM(basisStates, groundState, reducingIndices; reducingConfigs=reducingConfigs)
    eigenvalues = eigvals(Hermitian(reducedDMatrix))
    eigenvalues[eigenvalues .< 1e-16] .= 0
    @assert all(x -> x ≥ 0, eigenvalues)
    return -sum(eigenvalues .* log.(replace(eigenvalues, 0 => 1e-10)))
end


function mutInfo(basisStates::Vector{BitArray}, 
        groundState::Vector{Float64}, 
        reducingIndices::Tuple{Vector{Int64}, Vector{Int64}};
        reducingConfigs::Tuple{Vector{BitVector}, Vector{BitVector}}=(BitVector[], BitVector[])
    )
    combinedConfigs = vec([[c1; c2] for (c1, c2) in Iterators.product(reducingConfigs...)])
    SEE_A = vnEntropy(basisStates, groundState, reducingIndices[1]; reducingConfigs=reducingConfigs[1])
    SEE_B = vnEntropy(basisStates, groundState, reducingIndices[2]; reducingConfigs=reducingConfigs[2])
    SEE_AB = vnEntropy(basisStates, groundState, vcat(reducingIndices...); reducingConfigs=combinedConfigs)
    println((SEE_A, SEE_B, SEE_AB))
    return SEE_A + SEE_B - SEE_AB
end
