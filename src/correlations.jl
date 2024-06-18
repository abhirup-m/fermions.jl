function gstateCorrelation(gstates::Vector{Dict{BitVector, Float64}}, correlationOperator::Dict{Tuple{String,Vector{Int64}},Float64})
    correlationResult = 0
    for state in gstates
        @assert sum(collect(values(state)) .^ 2) ≈ 1
        opOnState = applyOperatorOnState(state, correlationOperator)
        stateCorr = 0
        for (bstate, coeff) in opOnState
            stateCorr += bstate ∈ keys(state) ? coeff * state[bstate] : 0
        end
        correlationResult += abs(stateCorr) / length(gstates)
    end
    return correlationResult
end


function vnEntropy(groundState::Dict{BitVector, Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[])
    reducedDMatrix = reducedDM(groundState, reducingIndices; reducingConfigs=reducingConfigs)
    eigenvalues = eigvals(Hermitian(reducedDMatrix))
    eigenvalues[eigenvalues .< 1e-16] .= 0
    @assert all(x -> x ≥ 0, eigenvalues)
    return -sum(eigenvalues .* log.(replace(eigenvalues, 0 => 1e-10)))
end


function mutInfo(
        groundState::Dict{BitVector, Float64}, 
        reducingIndices::Tuple{Vector{Int64}, Vector{Int64}};
        reducingConfigs::Tuple{Vector{BitVector}, Vector{BitVector}}=(BitVector[], BitVector[])
    )
    combinedConfigs = vec([[c1; c2] for (c1, c2) in Iterators.product(reducingConfigs...)])
    SEE_A = vnEntropy(groundState, reducingIndices[1]; reducingConfigs=reducingConfigs[1])
    SEE_B = vnEntropy(groundState, reducingIndices[2]; reducingConfigs=reducingConfigs[2])
    SEE_AB = vnEntropy(groundState, vcat(reducingIndices...); reducingConfigs=combinedConfigs)
    println((SEE_A, SEE_B, SEE_AB))
    return SEE_A + SEE_B - SEE_AB
end
