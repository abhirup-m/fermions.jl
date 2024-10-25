"""
    GenCorrelation(state, operator)

Calculate the expectation value ⟨state|operator|state⟩.

# Examples
```jldoctest
julia> state = Dict(Bool.([1, 0]) => 0.5, Bool.([0, 1]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 0] => 0.5
  [0, 1] => -0.5

julia> operator = [("+-", [1, 2], 1.0)]
1-element Vector{Tuple{String, Vector{Int64}, Float64}}:
 ("+-", [1, 2], 1.0)

julia> GenCorrelation(state, operator)
-0.5
```
"""
function GenCorrelation(
        state::Dict{BitVector,Float64},
        operator::Vector{Tuple{String,Vector{Int64},Float64}}
    )
    # state is of the form {|1>: c_1, |2>: c_2, ... |n>: c_n}.
    intermStates = fetch.([Threads.@spawn ApplyOperator([term], state) for term in operator])
    correlation = sum(fetch.([Threads.@spawn StateOverlap(state, intermState) for intermState in intermStates]))
    return correlation / sum(values(state) .^ 2)
end
export GenCorrelation


function GenCorrelation(
        state::Vector{Float64},
        operator::Matrix{Float64}
    )
    @assert length(state) == size(operator)[1]
    state ./= norm(state)
    return state' * operator * state
end
export GenCorrelation


"""
    reducedDM(state, reducingIndices)

Reduce the density matrix |state⟩⟨state| by tracing over the 
indices _not_ specified in reducingIndices.

# Examples
```jldoctest
julia> state = Dict(Bool.([1, 0]) => 0.5, Bool.([0, 1]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 0] => 0.5
  [0, 1] => -0.5

julia> reducedDM(state, [1])
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.5
```
"""
function ReducedDM(
        state::Dict{BitVector,Float64},
        reducingIndices::Vector{Int64};
        reducingConfigs::Vector{BitVector}=BitVector[]
    )

    # indices of the degrees of freedom that will not be summed over.
    nonReducingIndices = setdiff(1:length(collect(keys(state))[1]), reducingIndices)

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
    for (psi_im, c_im) in state

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

    return reducedDMatrix ./ sum(diag(reducedDMatrix))
end
export ReducedDM


"""
    VonNEntropy(state, reducingIndices)

Calculate entanglement entropy of the subsystem defined by
`reducingIndices`.

# Examples
```jldoctest
julia> state = Dict(Bool.([1, 0]) => 0.5, Bool.([0, 1]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 0] => 0.5
  [0, 1] => -0.5

julia> VonNEntropy(state, [1])
0.6931471805599453
```
"""
function VonNEntropy(
        state::Dict{BitVector,Float64},
        reducingIndices::Vector{Int64};
        reducingConfigs::Vector{BitVector}=BitVector[],
        tolerance=1e-10, schmidtGap=false
    )
    reducedDMatrix = ReducedDM(state, reducingIndices; reducingConfigs=reducingConfigs)
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
export VonNEntropy


"""
    MutInfo(state, reducingIndices)

Calculate mutual information between the subsystems defined by
the tuple `reducingIndices`.

# Examples
```jldoctest
julia> state = Dict(Bool.([1, 0, 1]) => 0.5, Bool.([0, 1, 0]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [0, 1, 0] => -0.5
  [1, 0, 1] => 0.5

julia> MutInfo(state, ([1], [2]))
0.6931471805599453
```
"""
function MutInfo(
    state::Dict{BitVector,Float64},
    reducingIndices::Tuple{Vector{Int64},Vector{Int64}};
    reducingConfigs::Tuple{Vector{BitVector},Vector{BitVector}}=(BitVector[], BitVector[])
)
    combinedConfigs = vec([[c1; c2] for (c1, c2) in Iterators.product(reducingConfigs...)])
    SEE_A = Threads.@spawn VonNEntropy(state, reducingIndices[1]; reducingConfigs=reducingConfigs[1])
    SEE_B = Threads.@spawn VonNEntropy(state, reducingIndices[2]; reducingConfigs=reducingConfigs[2])
    SEE_AB = Threads.@spawn VonNEntropy(state, vcat(reducingIndices...); reducingConfigs=combinedConfigs)
    return sum([1, 1, -1] .* fetch.([SEE_A, SEE_B, SEE_AB]))
end
export MutInfo


"""
    TripartiteInfo(state, reducingIndices)

Calculate tripartite information between the subsystems defined by
the 3-tuple `reducingIndices`.

# Examples
```jldoctest
julia> state = Dict(Bool.([1, 0, 1, 0]) => 0.5, Bool.([0, 1, 0, 1]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 0, 1, 0] => 0.5
  [0, 1, 0, 1] => -0.5

julia> TripartiteInfo(state, ([1], [2], [3]))
0.6931471805599453
```
"""
function TripartiteInfo(
    groundState::Dict{BitVector,Float64},
    reducingIndices::NTuple{3,Vector{Int64}};
    reducingConfigs::NTuple{3,Vector{BitVector}}=(BitVector[], BitVector[], BitVector[])
)
    BC_indices = vcat(reducingIndices[2], reducingIndices[3])
    BC_configs = vcat(reducingConfigs[2], reducingConfigs[3])
    I2_A_B = Threads.@spawn MutInfo(groundState, reducingIndices[[1, 2]]; reducingConfigs=reducingConfigs[[1, 2]])
    I2_A_C = Threads.@spawn MutInfo(groundState, reducingIndices[[1, 3]]; reducingConfigs=reducingConfigs[[1, 3]])
    I2_A_BC = Threads.@spawn MutInfo(groundState, (reducingIndices[1], BC_indices); reducingConfigs=(reducingConfigs[1], BC_configs))
    return sum([1, 1, -1] .* fetch.([I2_A_B, I2_A_C, I2_A_BC]))
end
export TripartiteInfo


"""
    ThermalAverage(eigenStates, eigenVals, operator, invTemp)

Calculate the canonical ensemble average of `operator` at the given
inverse temperature for the given spectrum.

# Examples
```jldoctest
julia> basis = BasisStates(2)
4-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 0] => 1.0)
 Dict([0, 1] => 1.0)
 Dict([1, 0] => 1.0)
 Dict([1, 1] => 1.0)

julia> eigenStates = BasisStates(2)
4-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 0] => 1.0)
 Dict([0, 1] => 1.0)
 Dict([1, 0] => 1.0)
 Dict([1, 1] => 1.0)

julia> eigenVals = rand(4)
4-element Vector{Float64}:
 0.018726877619070992
 0.43796674226812815
 0.7905965730299785
 0.97571966188274

julia> operator = [("n", [1], 1.0)]
1-element Vector{Tuple{String, Vector{Int64}, Float64}}:
 ("n", [1], 1.0)

julia> ThermalAverage(basis, eigenVals, operator, 1.0)
0.33797199716422116
```
"""
function ThermalAverage(
    eigenStates::Vector{Dict{BitVector,Float64}},
    eigenVals::Vector{Float64},
    operator::Vector{Tuple{String,Vector{Int64},Float64}},
    invTemp::Float64,
)
    return sum(fetch.([Threads.@spawn exp(-invTemp * energy) * GenCorrelation(state, operator) for (state, energy) in zip(eigenStates, eigenVals)])) / sum(exp.(-invTemp .* eigenVals))
end
export ThermalAverage


function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Vector{Float64}},
    probes::Dict{String,Matrix{Float64}},
    freqValues::Vector{Float64},
    standDev::Float64;
    )
    @assert length(eigVals) == length(eigVecs)

    groundState = eigVecs[sortperm(eigVals)[1]]
    energyGs = minimum(eigVals)
    specFunc = 0 .* freqValues
    spectralWeights = [(0.,0.) for _ in eigVecs]
    @showprogress Threads.@threads for index in eachindex(eigVecs)
        excitedState = eigVecs[index]
        spectralWeights[index] = ((groundState' * probes["destroy"] * excitedState) * (excitedState' * probes["create"] * groundState), 
                                  (excitedState' * probes["destroy"] * groundState) * (groundState' * probes["create"] * excitedState)
                                 )
    end
    for index in eachindex(eigVals)
        specFunc .+= spectralWeights[index][1] * standDev ./ ((freqValues .+ energyGs .- eigVals[index]) .^ 2 .+ standDev^2)
        specFunc .+= spectralWeights[index][2] * standDev ./ ((freqValues .- energyGs .+ eigVals[index]) .^ 2 .+ standDev^2)
    end
    return specFunc ./ sum(specFunc .* (maximum(freqValues) - minimum(freqValues[1])) / (length(freqValues) - 1))
end
export SpecFunc


"""
    SpecFunc(eigVals, eigVecs, probe, probeDiag, freqValues, standDev)

Calculate the spectral function for the excitations defined by probe
and probeDiag.

# Examples
```jldoctest
julia> eigenVals = [0., 1., 1.];

julia> eigenStates = Dict{BitVector, Float64}[Dict([1, 0] => 1.0, [0, 1] => 1.0), Dict([1, 1] => 1.0), Dict([0, 0] => 1.0)]
3-element Vector{Dict{BitVector, Float64}}:
 Dict([1, 0] => 1.0, [0, 1] => 1.0)
 Dict([1, 1] => 1.0)
 Dict([0, 0] => 1.0)

julia> probe = [("-", [1], 1.0)];

julia> probeDag = [("+", [1], 1.0)];

julia> freqValues = collect(range(-2, stop=2, length=10));

julia> SpecFunc(eigenVals, eigenStates, probe, probeDag, freqValues, 1e-2)
10-element Vector{Float64}:
 0.011110098865559272
 0.03392067328129922
 0.8057354340607891
 0.09351894323911442
 0.023221646866030225
 0.023221646866030225
 0.09351894323911442
 0.8057354340607891
 0.03392067328129922
 0.011110098865559272
```
"""
function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probes::Dict{String,Vector{Tuple{String,Vector{Int64},Float64}}},
    freqValues::Vector{Float64},
    basisStates::Vector{Dict{BitVector,Float64}},
    standDev::Float64;
)

    eigenStates = [ExpandIntoBasis(vector, basisStates) for vector in eigVecs]

    probes = Dict{String,Matrix{Float64}}(name => OperatorMatrix(basisStates, probe) for (name,probe) in probes)

    return SpecFunc(eigVals, eigenStates, probes, freqValues, standDev)
end
export SpecFunc


"""
    SpecFunc(eigVals, eigVecs, probe, probeDiag, freqValues, standDev, symmetries)

Extends SpecFunc() by making use of symmetries.

# Examples
```jldoctest
julia> SpecFunc(eigenVals, eigenStates, probe, probeDag, freqValues, 1e-2, ['N'])
```
"""
function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probes::Dict{String,Vector{Tuple{String,Vector{Int64},Float64}}},
    freqValues::Vector{Float64},
    basisStates::Vector{Dict{BitVector,Float64}},
    standDev::Float64,
    symmetries::Vector{Char};
)

    groundState = eigVecs[sortperm(eigVals)[1]]

    # calculate sector of excited states
    excitedSectorHole = GetSector(ApplyOperator(probes["destroy"], groundState), symmetries)
    excitedSectorParticle = GetSector(ApplyOperator(probes["create"], groundState), symmetries)

    classifiedSpectrum, classifiedEnergies = ClassifyBasis(eigVecs, symmetries; energies=eigVals)
    minimalEigVecs = Dict{BitVector, Float64}[[groundState]; classifiedSpectrum[excitedSectorHole]; classifiedSpectrum[excitedSectorParticle];]
    minimalEigVals = Float64[[minimum(eigVals)]; classifiedEnergies[excitedSectorHole]; classifiedEnergies[excitedSectorParticle];]

    return SpecFunc(minimalEigVals, minimalEigVecs, probes, freqValues, basisStates, standDev)
end
export SpecFunc


function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probes::Dict{String,Vector{Tuple{String,Vector{Int64},Float64}}},
    freqValues::Vector{Float64},
    basisStates::Vector{Dict{BitVector,Float64}},
    standDev::Float64,
    symmetries::Vector{Char},
    groundStateSector::Union{Tuple{Int64}, Tuple{Int64,Int64}},
)

    @assert length(eigVals) == length(eigVecs)
    classifiedSpectrum, classifiedEnergies = ClassifyBasis(eigVecs, symmetries; energies=eigVals)

    groundState = classifiedSpectrum[groundStateSector][sortperm(classifiedEnergies[groundStateSector])[1]]

    # calculate sector of excited states
    excitedSectorHole = GetSector(ApplyOperator(probes["destroy"], groundState), symmetries)
    excitedSectorParticle = GetSector(ApplyOperator(probes["create"], groundState), symmetries)

    minimalEigVecs = Dict{BitVector, Float64}[[groundState]; classifiedSpectrum[excitedSectorHole]; classifiedSpectrum[excitedSectorParticle];]
    minimalEigVals = Float64[[minimum(classifiedEnergies[groundStateSector])]; classifiedEnergies[excitedSectorHole]; classifiedEnergies[excitedSectorParticle];]
    @assert minimalEigVals[1] == minimum(minimalEigVals)

    return SpecFunc(minimalEigVals, minimalEigVecs, probes, freqValues, basisStates, standDev)
end
export SpecFunc


"""
    SpecFunc(eigVals, eigVecMatrix, probe, probeDiag, freqValues, standDev)

Extends SpecFunc() by directly accepting matrices instead of abstract states and operators.
Useful when the matrix representations are known but the basis is very complicated.
"""
function SpecFunc(
    eigVals::Vector{Float64},
    eigVecMatrix::Matrix{Float64},
    probes::Dict{String,Matrix{Float64}},
    freqValues::Vector{Float64},
    standDev::Float64,
)

    eigVecs = [collect(vec) for vec in eachcol(eigVecMatrix)]
    return SpecFunc(eigVals, eigVecs, probes, freqValues, standDev)
end
export SpecFunc
