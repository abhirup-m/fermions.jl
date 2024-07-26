"""Calculate static correlation function ⟨ψ|O|ψ⟩ of the provided operator."""
function GenCorrelation(state::Dict{BitVector,Float64}, operator::Vector{Tuple{String,Vector{Int64},Float64}})
    # Gstate vector is of the form {|1>: c_1, |2>: c_2, ... |n>: c_n}.
    intermStates = fetch.([Threads.@spawn ApplyOperator([term], state) for term in operator])
    correlation = sum(fetch.([Threads.@spawn StateOverlap(state, intermState) for intermState in intermStates]))
    return correlation / sum(values(state) .^ 2)
end


"""Given a state |Ψ>, calculates the reduced density matrix ρ_A = ∑_i ⟨i|ρ|i⟩ for the subsystem A where 
i sums over all configurations of the subsystem A. The subsystem A is specified by the argument
`reducingIndices` which specifies the indices in the representation that go into forming A. Also
accepts an optional argument `reducingConfigs` that specifies a subset of states of A to sum over.
"""
function reducedDM(groundState::Dict{BitVector,Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[])

    # indices of the degrees of freedom that will not be summed over.
    nonReducingIndices = setdiff(1:length(collect(keys(groundState))[1]), reducingIndices)

    # get all the configurations of the subsystem A in order to perform the partial trace.
    if length(reducingConfigs) == 0
        reducingConfigs = [collect(config) for config in Iterators.product(fill([1, 0], length(reducingIndices))...)]
    end
    reducingConfigsMaps = Dict(config => i for (i, config) in enumerate(reducingConfigs))

    # reducingCoeffs = Dict{BitVector,Dict{BitVector,Float64}}(BitVector(config) => Dict()
    # for config in reducingConfigs)
    #
    reducingCoeffs = Dict{BitVector,Vector{Float64}}()
    reducingIndexSets = Dict{BitVector,Vector{Int64}}()
    reducedDMatrix = zeros(length(reducingConfigs), length(reducingConfigs))

    # |ψ⟩ can be written as ∑ c_{im} |i>|m>, where i is a state in the subspace A,
    # while m is a state in the subspace A-complement.
    for (psi_im, c_im) in groundState

        # get |i>
        psi_i = psi_im[reducingIndices]
        index_i = reducingConfigsMaps[psi_i]

        # get |m>
        psi_m = psi_im[nonReducingIndices]

        if psi_i ∉ reducingConfigs
            continue
        end

        if psi_m ∉ keys(reducingCoeffs)
            reducingCoeffs[psi_m] = [c_im]
            reducingIndexSets[psi_m] = [index_i]
        else
            reducedDMatrix[index_i, reducingIndexSets[psi_m]] .+= c_im .* reducingCoeffs[psi_m]
            reducedDMatrix[reducingIndexSets[psi_m], index_i] .+= c_im .* reducingCoeffs[psi_m]
            push!(reducingCoeffs[psi_m], c_im)
            push!(reducingIndexSets[psi_m], reducingConfigsMaps[psi_i])
        end
        reducedDMatrix[index_i, index_i] += c_im^2


    end

    return reducedDMatrix
end


function vnEntropy(groundState::Dict{BitVector,Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[], tolerance=1e-10, schmidtGap=false)
    reducedDMatrix = reducedDM(groundState, reducingIndices; reducingConfigs=reducingConfigs)
    eigenvalues = eigvals(0.5 * (reducedDMatrix + reducedDMatrix'))
    eigenvalues[eigenvalues.<tolerance] .= 0
    eigenvalues ./= sum(eigenvalues)
    @assert all(x -> x ≥ 0, eigenvalues)
    @assert isapprox(sum(eigenvalues), 1; atol=tolerance)
    if !schmidtGap
        return -sum(eigenvalues .* log.(replace(eigenvalues, 0 => 1e-10)))
    else
        largestEigvals = sort(eigenvalues)
        return -sum(eigenvalues .* log.(replace(eigenvalues, 0 => 1e-10))), largestEigvals[end] - largestEigvals[end-1]
    end
end


function mutInfo(
    groundState::Dict{BitVector,Float64},
    reducingIndices::Tuple{Vector{Int64},Vector{Int64}};
    reducingConfigs::Tuple{Vector{BitVector},Vector{BitVector}}=(BitVector[], BitVector[])
)
    combinedConfigs = vec([[c1; c2] for (c1, c2) in Iterators.product(reducingConfigs...)])
    SEE_A = Threads.@spawn vnEntropy(groundState, reducingIndices[1]; reducingConfigs=reducingConfigs[1])
    SEE_B = Threads.@spawn vnEntropy(groundState, reducingIndices[2]; reducingConfigs=reducingConfigs[2])
    SEE_AB = Threads.@spawn vnEntropy(groundState, vcat(reducingIndices...); reducingConfigs=combinedConfigs)
    return sum([1, 1, -1] .* fetch.([SEE_A, SEE_B, SEE_AB]))
end


function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probe::Vector{Tuple{String,Vector{Int64},Float64}},
    probeDag::Vector{Tuple{String,Vector{Int64},Float64}},
    freqArray::Vector{Float64},
    broadening::Float64
)

    energyGs, groundState = eigVals[1], eigVecs[1]
    gsBasisState = collect(keys(groundState))[1]
    excitedBasisState = collect(keys(ApplyOperator(probe, groundState)))[1]
    excitedBasisStateDag = collect(keys(ApplyOperator(probeDag, groundState)))[1]

    excitedSectors = [
                      (sum(excitedBasisState), sum(excitedBasisState[1:2:end]) - sum(excitedBasisState[2:2:end])),
                      (sum(excitedBasisStateDag), sum(excitedBasisStateDag[1:2:end]) - sum(excitedBasisStateDag[2:2:end]))
                     ]

    excitedVecs = Dict(excitedSectors .=> [Dict{BitVector, Float64}[], Dict{BitVector, Float64}[]])
    excitedVals = Dict(excitedSectors .=> [Float64[], Float64[]])
    for (energy, state) in zip(eigVals, eigVecs)
        basisState = collect(keys(state))[1]
        sector = (sum(basisState), sum(basisState[1:2:end]) - sum(basisState[2:2:end]))
        if sector in excitedSectors
            push!(excitedVecs[sector], state)
            push!(excitedVals[sector], energy)
        end
    end

    specFuncArray = SpecFunc((energyGs, groundState), excitedVals, excitedVecs, probe, probeDag, freqArray, broadening)
    return specFuncArray
end


function SpecFunc(
    groundSE::Tuple{Float64,Dict{BitVector,Float64}},
    eigVals::Dict{Tuple{Int64,Int64},Vector{Float64}},
    eigVecs::Dict{Tuple{Int64,Int64},Vector{Dict{BitVector,Float64}}},
    probe::Vector{Tuple{String,Vector{Int64},Float64}},
    probeDag::Vector{Tuple{String,Vector{Int64},Float64}},
    freqArray::Vector{Float64},
    broadening::Float64
)

    energyGs, groundState = groundSE

    # calculate c_ν |GS>
    #=display(sort([v for (k, v) in groundState if k[1] == 1], rev=true))=#
    excitedState = ApplyOperator(probe, groundState)
    exampleState = collect(keys(excitedState))[1]
    excitedSector = (sum(exampleState), sum(exampleState[1:2:end] - exampleState[2:2:end]))

    # calculate c^†_ν |GS>
    excitedStateDag = ApplyOperator(probeDag, groundState)
    exampleState = collect(keys(excitedStateDag))[1]
    excitedSectorDag = (sum(exampleState), sum(exampleState[1:2:end] - exampleState[2:2:end]))

    # create array of frequency points and spectral function
    specFuncArrayP = [0 .* freqArray for _ in eigVals[excitedSector]]
    specFuncArrayH = [0 .* freqArray for _ in eigVals[excitedSectorDag]]

    # loop over eigenstates of the excited symmetry sectors,
    # calculated overlaps and multiply the appropriate denominators.
    @sync begin
        @async Threads.@threads for index in eachindex(eigVals[excitedSector])
            particleWeight = StateOverlap(eigVecs[excitedSector][index], excitedState)^2
            specFuncArrayP[index] = particleWeight * broadening ./ ((freqArray .- energyGs .+ eigVals[excitedSector][index]) .^ 2 .+ broadening^2)
        end
        @async Threads.@threads for index in eachindex(eigVals[excitedSectorDag])
            holeWeight = StateOverlap(eigVecs[excitedSectorDag][index], excitedStateDag)^2
            specFuncArrayH[index] = holeWeight * broadening ./ ((freqArray .+ energyGs .- eigVals[excitedSectorDag][index]) .^ 2 .+ broadening^2)
        end
    end
    specFuncArray = sum(specFuncArrayP) .+ sum(specFuncArrayH)
    return specFuncArray
end


function tripartiteInfo(
    groundState::Dict{BitVector,Float64},
    reducingIndices::NTuple{3,Vector{Int64}};
    reducingConfigs::NTuple{3,Vector{BitVector}}=(BitVector[], BitVector[], BitVector[])
)
    BC_indices = vcat(reducingIndices[2], reducingIndices[3])
    BC_configs = vcat(reducingConfigs[2], reducingConfigs[3])
    I2_A_B = Threads.@spawn mutInfo(groundState, reducingIndices[[1, 2]]; reducingConfigs=reducingConfigs[[1, 2]])
    I2_A_C = Threads.@spawn mutInfo(groundState, reducingIndices[[1, 3]]; reducingConfigs=reducingConfigs[[1, 3]])
    I2_A_BC = Threads.@spawn mutInfo(groundState, (reducingIndices[1], BC_indices); reducingConfigs=(reducingConfigs[1], BC_configs))
    return sum([1, 1, -1] .* fetch.([I2_A_B, I2_A_C, I2_A_BC]))
end


function ThermalAverage(
    eigenStates::Vector{Dict{BitVector,Float64}},
    eigenVals::Vector{Float64},
    operator::Vector{Tuple{String,Vector{Int64},Float64}},
    invTemp::Float64,
)
    return sum(fetch.([Threads.@spawn exp(-invTemp * energy) * GenCorrelation(state, operator) for (state, energy) in zip(eigenStates, eigenVals)])) / sum(exp.(-invTemp .* eigenVals))
end
