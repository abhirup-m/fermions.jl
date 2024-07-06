"""Calculate static correlation function ⟨ψ|O|ψ⟩ of the provided operator."""
function gstateCorrelation(gstate::Dict{BitVector, Float64}, correlationOperator::Dict{Tuple{String,Vector{Int64}},Float64})
    # Gstate vector is of the form {|1>: c_1, |2>: c_2, ... |n>: c_n}.
    # loop over the pairs (|m>, c_m)
    
    # check that state is normalised: ∑ c_m^2 = 1
    @assert isapprox(sum(values(gstate) .^ 2), 1, atol = 1e-6)

    # check the action of the operator O on the ground state
    opOnState = applyOperatorOnState(gstate, correlationOperator)

    # calculate the overlap of the new state with the ground state
    expecValue = sum([gstate[key] * opOnState[key] for key in intersect(keys(opOnState), keys(gstate))])
    return expecValue
end


function reducedDM(groundState::Dict{BitVector, Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[])
    if length(reducingConfigs) == 0
        reducingConfigs = [collect(config) for config in Iterators.product(fill([1, 0], length(reducingIndices))...)]
    end

    nonReducingIndices = setdiff(1:length(collect(keys(groundState))[1]), reducingIndices)
    reducingCoeffs = Dict{BitVector, Dict{BitVector, Float64}}(BitVector(config) => Dict()
                                                                     for config in reducingConfigs)
    for (state, coeff) in groundState
        reducingConfig = state[reducingIndices]
        nonReducingConfig = state[nonReducingIndices]
        if reducingConfig ∉ reducingConfigs
            continue
        end
        reducingCoeffs[reducingConfig][nonReducingConfig] = coeff
    end
    reducedDMatrix = zeros(length(reducingConfigs), length(reducingConfigs))
    Threads.@threads for ((i1, c1), (i2, c2)) in collect(Iterators.product(enumerate(reducingConfigs), enumerate(reducingConfigs)))
        reducedDMatrix[i1, i2] = sum(values(merge(*, reducingCoeffs[collect(c1)], reducingCoeffs[collect(c2)])))
    end
    return reducedDMatrix
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
        eigVals::Tuple{Dict{Tuple{Int64, Int64}, Vector{Float64}}, Dict{Tuple{Int64, Int64}, Vector{Float64}}}, 
        eigVecs::Tuple{Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}}, Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}}},
        probe::Dict{Tuple{String,Vector{Int64}},Float64},
        probeDag::Dict{Tuple{String,Vector{Int64}},Float64},
        freqArray::Vector{Float64},
        broadening::Float64
    )

    # calculate c_ν |GS>
    excitedState = applyOperatorOnState(groundState, probe)

    # calculate c^†_ν |GS>
    excitedStateDag = applyOperatorOnState(groundState, probeDag)

    # calculate the quantum numbers of the symmetry sector in which
    # the above excited states reside. We only need to the check the
    # eigenstates in these symmetry sectors to calculate the overlaps.
    excCompState = collect(keys(excitedState))[1]
    excitedSector = (sum(excCompState), sum(excCompState[1:2:end]) - sum(excCompState[2:2:end]))
    excCompStateDag = collect(keys(excitedStateDag))[1]
    excitedSectorDag = (sum(excCompStateDag), sum(excCompStateDag[1:2:end]) - sum(excCompStateDag[2:2:end]))

    # create array of frequency points and spectral function
    specFuncArray = 0 .* freqArray

    # loop over eigenstates of the excited symmetry sectors,
    # calculated overlaps and multiply the appropriate denominators.
    for (i, stateDict) in collect(enumerate(eigVecs[1][excitedSector]))
        for key in intersect(keys(excitedState), keys(stateDict))
            overlap = abs(stateDict[key] * excitedState[key])^2
            specFuncArray .+= overlap * broadening ./ ((freqArray .+ energyGs .- eigVals[1][excitedSector][i]) .^ 2 .+ broadening ^ 2)
        end
    end
    for (i, stateDict) in collect(enumerate(eigVecs[2][excitedSectorDag]))
        for key in intersect(keys(excitedStateDag), keys(stateDict))
            overlapDag = abs(stateDict[key] * excitedStateDag[key])^2
            specFuncArray .+= overlapDag * broadening ./ ((freqArray .- energyGs .+ eigVals[2][excitedSectorDag][i]) .^ 2 .+ broadening ^ 2)
        end
    end
    specFuncArray ./= sum(specFuncArray) * abs(freqArray[2] - freqArray[1])
    return specFuncArray
end
