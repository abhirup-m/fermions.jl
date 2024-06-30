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


function specFunc(
        groundState::Dict{BitVector,Float64},
        energyGs::Float64,
        eigVals::Dict{Tuple{Int64, Int64}, Vector{Float64}}, 
        eigVecs::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}},
        probe::Dict{Tuple{String,Vector{Int64}},Float64},
        probeDag::Dict{Tuple{String,Vector{Int64}},Float64},
        freqPoints::Tuple{Float64, Int64},
        broadening::Float64
    )

    # calculate c_ν |GS>
    excitedState = applyOperatorOnState(groundState, probe)

    # calculate c^†_ν |GS>
    excitedStateDag = applyOperatorOnState(groundState, probeDag)

    # calculate the quantum numbers of the symmetry sector in which
    # the above excited states reside. We only need to the check the
    # eigenstates in these symmetry sectors to calculate the overlaps.
    excitedSector = (sum(collect(keys(excitedState))[1]), 
                     sum(collect(keys(excitedState))[1:2:end]) - sum(collect(keys(excitedState))[2:2:end])
                    )
    excitedSectorDag = (sum(collect(keys(excitedStateDag))[1]), 
                     sum(collect(keys(excitedStateDag))[1:2:end]) - sum(collect(keys(excitedStateDag))[2:2:end])
                    )

    # create array of frequency points and spectral function
    freqArray = range(-abs(freqPoints[1]), stop=abs(freqPoints[1]), length=freqPoints[2])
    specFuncArray = 0 .* freqArray

    # loop over eigenstates of the excited symmetry sectors,
    # calculated overlaps and multiply the appropriate denominators.
    for (i, stateDict) in collect(enumerate(eigVecs[excitedSector]))
        overlap = sum([stateDict[key] * excitedState[key] for key in intersect(keys(excitedState), keys(stateDict))])
        specFuncArray .+= abs(overlap)^2 * broadening ./ ((freqArray + energyGs .- eigVals[key][i]) .^ 2 .+ broadening ^ 2)
    end
    for (i, stateDict) in collect(enumerate(eigVecs[excitedSectorDag]))
        overlapDag = sum([stateDict[key] * excitedStateDag[key] for key in intersect(keys(excitedStateDag), keys(stateDict))])
        specFuncArray .+= abs(overlapDag)^2 * broadening ./ ((freqArray - energyGs .+ eigVals[excitedSectorDag][i]) .^ 2 .+ broadening ^ 2)
    end
    return specFuncArray
end
