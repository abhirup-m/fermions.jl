"""Calculate static correlation function ⟨ψ|O|ψ⟩ of the provided operator."""
function FastCorrelation(state::Dict{BitVector,Float64}, operator::Vector{Tuple{String,Vector{Int64},Float64}})
    # Gstate vector is of the form {|1>: c_1, |2>: c_2, ... |n>: c_n}.
    # loop over the pairs (|m>, c_m)

    # check that state is normalised: ∑ c_m^2 = 1
    @assert isapprox(sum(values(state) .^ 2), 1, atol=1e-6)


    correlation = sum(fetch.([Threads.@spawn StateOverlap(state, ApplyOperator([term], state))
                              for term in operator]))
    return correlation
end


"""Given a state |Ψ>, calculates the reduced density matrix ρ_A = ∑_i ⟨i|ρ|i⟩ for the subsystem A where 
i sums over all configurations of the subsystem A. The subsystem A is specified by the argument
`reducingIndices` which specifies the indices in the representation that go into forming A. Also
accepts an optional argument `reducingConfigs` that specifies a subset of states of A to sum over.
"""
function reducedDM(groundState::Dict{BitVector,Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[])

    # get all the configurations of the subsystem A in order to perform the partial trace.
    if length(reducingConfigs) == 0
        reducingConfigs = [collect(config) for config in Iterators.product(fill([1, 0], length(reducingIndices))...)]
    end

    # indices of the degrees of freedom that will not be summed over.
    nonReducingIndices = setdiff(1:length(collect(keys(groundState))[1]), reducingIndices)
    reducingCoeffs = Dict{BitVector,Dict{BitVector,Float64}}(BitVector(config) => Dict()
                                                             for config in reducingConfigs)

    # |ψ⟩ can be written as ∑ c_{im} |i>|m>, where i is a state in the subspace A,
    # while m is a state in the subspace A-complement.
    for (state, coeff) in groundState

        # get |i>
        reducingConfig = state[reducingIndices]

        # get |m>
        nonReducingConfig = state[nonReducingIndices]
        if reducingConfig ∉ reducingConfigs
            continue
        end

        # store c_im by classifying according to i and m
        reducingCoeffs[reducingConfig][nonReducingConfig] = coeff
    end

    reducedDMatrix = zeros(length(reducingConfigs), length(reducingConfigs))
    Threads.@threads for ((i1, c1), (i2, c2)) in collect(Iterators.product(enumerate(reducingConfigs), enumerate(reducingConfigs)))
        commonKeys = intersect(keys(reducingCoeffs[collect(c1)]), keys(reducingCoeffs[collect(c2)]))
        reducedDMatrix[i1, i2] = sum([reducingCoeffs[collect(c1)][key] * reducingCoeffs[collect(c2)][key] for key in commonKeys])
    end
    return reducedDMatrix
end


function vnEntropy(groundState::Dict{BitVector,Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[], tolerance=1e-10)
    reducedDMatrix = reducedDM(groundState, reducingIndices; reducingConfigs=reducingConfigs)

    eigenvalues = eigvals(reducedDMatrix)

    eigenvalues[eigenvalues.<tolerance] .= 0
    eigenvalues ./= sum(eigenvalues)
    @assert all(x -> x ≥ 0, eigenvalues)
    @assert isapprox(sum(eigenvalues), 1; atol=tolerance)
    return -sum(eigenvalues .* log.(replace(eigenvalues, 0 => 1e-10)))
end


function mutInfo(
    groundState::Dict{BitVector,Float64},
    reducingIndices::Tuple{Vector{Int64},Vector{Int64}};
    reducingConfigs::Tuple{Vector{BitVector},Vector{BitVector}}=(BitVector[], BitVector[])
)
    combinedConfigs = vec([[c1; c2] for (c1, c2) in Iterators.product(reducingConfigs...)])
    SEE_A = vnEntropy(groundState, reducingIndices[1]; reducingConfigs=reducingConfigs[1])
    SEE_B = vnEntropy(groundState, reducingIndices[2]; reducingConfigs=reducingConfigs[2])
    SEE_AB = vnEntropy(groundState, vcat(reducingIndices...); reducingConfigs=combinedConfigs)
    return SEE_A + SEE_B - SEE_AB
end


function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probe::Tuple{String,Vector{Int64},Float64},
    probeDag::Tuple{String,Vector{Int64},Float64},
    freqArray::Vector{Float64},
    broadening::Float64
)

    groundState = eigVecs[1]
    energyGs = eigVals[1]

    # calculate c_ν |GS>
    excitedState = ApplyOperator([probe], groundState)

    # calculate c^†_ν |GS>
    excitedStateDag = ApplyOperator([probeDag], groundState)

    # create array of frequency points and spectral function
    specFuncArray = 0 .* freqArray

    # loop over eigenstates of the excited symmetry sectors,
    # calculated overlaps and multiply the appropriate denominators.
    for (eigvl, stateDict) in zip(eigVals, eigVecs)
        particleWeight = StateOverlap(stateDict, excitedState)^2
        holeWeight = StateOverlap(stateDict, excitedStateDag)^2

        specFuncArray .+= particleWeight * broadening ./ ((freqArray .- energyGs .+ eigvl) .^ 2 .+ broadening^2)
        specFuncArray .+= holeWeight * broadening ./ ((freqArray .+ energyGs .- eigvl) .^ 2 .+ broadening^2)
    end
    specFuncArray ./= sum(specFuncArray) * abs(freqArray[2] - freqArray[1])
    return specFuncArray
end
