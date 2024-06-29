"""Calculate static correlation function ⟨ψ|O|ψ⟩ of the provided operator."""
function gstateCorrelation(gstates::Dict{BitVector, Float64}, correlationOperator::Dict{Tuple{String,Vector{Int64}},Float64})
    # Gstate vector is of the form {|1>: c_1, |2>: c_2, ... |n>: c_n}.
    # loop over the pairs (|m>, c_m)
    
    # check that state is normalised: ∑ c_m^2 = 1
    @assert sum(values(state) .^ 2) ≈ 1

    # check the action of the operator O on the ground state
    opOnState = applyOperatorOnState(state, correlationOperator)

    # calculate the overlap of the new state with the ground state
    expecValue = sum([state[key] * opOnState[key] for key in intersect(keys(opOnState), keys(state))])
    return expecValue
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
    return SEE_A + SEE_B - SEE_AB
end
