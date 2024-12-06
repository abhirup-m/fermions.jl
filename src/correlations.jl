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


"""
    GenCorrelation(state, operator)

Calculates the expectation value <state|operator|state>, where state and operator
are accepted as vectors and matrices instead of dicts and bitvectors.
"""
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
        nonTracedSites::Vector{Int64};
        reducingConfigs::Vector{BitVector}=BitVector[]
    )
    # indices of the degrees of freedom that will not be summed over.
    tracedSites = setdiff(1:length(collect(keys(state))[1]), nonTracedSites)

    gatheringPermutation = vcat(nonTracedSites, tracedSites)
    permutedState = PermuteSites(deepcopy(state), gatheringPermutation)

    smallerBasis = [collect(keys(b))[1] for b in BasisStates(length(nonTracedSites))]
    reducedDMatrix = zeros(length(smallerBasis), length(smallerBasis))
    for (leftState, val1) in permutedState
        leftBasisState = Bool.(leftState[findall(∈(nonTracedSites), gatheringPermutation)])
        nonReducingConfig = leftState[findall(∈(tracedSites), gatheringPermutation)]
        for (rightState, val2) in permutedState
            if rightState[findall(∈(tracedSites), gatheringPermutation)] ≠ nonReducingConfig
                continue
            end
            rightBasisState = Bool.(rightState[findall(∈(nonTracedSites), gatheringPermutation)])
            leftIndex = findfirst(==(leftBasisState), smallerBasis)
            rightIndex = findfirst(==(rightBasisState), smallerBasis)
            reducedDMatrix[leftIndex, rightIndex] += val1 * val2'
        end
    end
    return reducedDMatrix ./ sum(diag(reducedDMatrix))
end
export ReducedDM


"""
    PartialTraceProjectors(nonTracedSites)

A reduced density matrix ρ_A for the subspace A can be calculate using the
relation ρ_A[i, j] = ⟨Ψ| |j⟩⟨i| |Ψ⟩, where |Ψ⟩ is of course the pure state
from which are calculating the reduced density matrix, and {|i⟩} are the
orthonormal basis states of subspace A. The object |j⟩⟨i| is therefore an
operator that, in general, transitions among various basis states. This 
function returns the complete set of transition operators {|j⟩⟨i|} for the
the given subspace A.
# Examples
```jldoctest
julia> operators = PartialTraceProjectors([2, 3]);

julia> println(join(operators, "\n"))
[("hhhh", [2, 3, 3, 2], 1.0)]
[("hh--", [2, 3, 3, 2], 1.0)]
[("hhh-", [2, 3, 3, 2], 1.0)]
[("hh-h", [2, 3, 3, 2], 1.0)]
[("++hh", [2, 3, 3, 2], 1.0)]
[("++--", [2, 3, 3, 2], 1.0)]
[("++h-", [2, 3, 3, 2], 1.0)]
[("++-h", [2, 3, 3, 2], 1.0)]
[("+hhh", [2, 3, 3, 2], 1.0)]
[("+h--", [2, 3, 3, 2], 1.0)]
[("+hh-", [2, 3, 3, 2], 1.0)]
[("+h-h", [2, 3, 3, 2], 1.0)]
[("h+hh", [2, 3, 3, 2], 1.0)]
[("h+--", [2, 3, 3, 2], 1.0)]
[("h+h-", [2, 3, 3, 2], 1.0)]
[("h+-h", [2, 3, 3, 2], 1.0)]
```
"""
function PartialTraceProjectors(
        nonTracedSites::Vector{Int64};
        nonTracedConfigs::Vector{BitVector}=BitVector[]
    )
    newBasis = mergewith(+, BasisStates(length(nonTracedSites))...) |> keys |> collect
    partialTraceProjectors = Array{Union{Nothing, Vector{Tuple{String, Vector{Int64}, Float64}}}}(nothing, length(newBasis), length(newBasis))
    for (leftIndex, rightIndex) in Iterators.product(eachindex(newBasis), eachindex(newBasis))
        leftBasisState = newBasis[leftIndex]
        rightBasisState = newBasis[rightIndex]
        preOperator = join(map(c -> c == 1 ? '+' : 'h', newBasis[rightIndex]))
        postOperator = join(map(c -> c == 1 ? '+' : 'h', newBasis[leftIndex]))
        postOperatorDag = Dagger(postOperator, nonTracedSites)
        partialTraceProjectors[leftIndex, rightIndex] = [(preOperator * postOperatorDag[1], vcat(nonTracedSites, postOperatorDag[2]), 1.)]
    end
    return partialTraceProjectors
end
export PartialTraceProjectors


"""
    ReducedDMProjectorBased(state, nonTracedSites)

Alternative method of calculating reduced density matrix. A reduced density 
matrix ρ_A for the subspace A can be calculate using the relation 
ρ_A[i, j] = ⟨Ψ| |j⟩⟨i| |Ψ⟩, where |Ψ⟩ is of course the pure state from which 
are calculating the reduced density matrix, and {|i⟩} are the orthonormal 
basis states of subspace A. The object |j⟩⟨i| is therefore an operator that,
in general, transitions among various basis states. 
"""
function ReducedDMProjectorBased(
        state::Dict{BitVector,Float64},
        nonTracedSites::Vector{Int64};
        nonTracedConfigs::Vector{BitVector}=BitVector[]
    )
    partialTraceProjectors = PartialTraceProjectors(nonTracedSites, nonTracedConfigs=nonTracedConfigs)
    reducedDMatrix = zeros(size(partialTraceProjectors)...)
    @assert size(partialTraceProjectors)[1] == size(partialTraceProjectors)[2]
    for leftIndex in 1:size(partialTraceProjectors)[1]
        for rightIndex in leftIndex:size(partialTraceProjectors)[2]
            reducedDMatrix[leftIndex, rightIndex] += GenCorrelation(state, partialTraceProjectors[leftIndex, rightIndex])
            reducedDMatrix[rightIndex, leftIndex] = reducedDMatrix[leftIndex, rightIndex]'
        end
    end
    return reducedDMatrix ./ sum(diag(reducedDMatrix))
end
export ReducedDMProjectorBased


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


"""
    SpecFunc(eigVals, eigVecs, probe, freqValues, standDev)

Fundamental implementation of spectral function calculation given a spectrum and excitation operators(defined by probe).
"""
function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Vector{Float64}},
    probes::Dict{String,Matrix{Float64}},
    freqValues::Vector{Float64},
    standDev::Union{Vector{Float64}, Float64};
    degenTol::Float64=0.,
    normalise::Bool=true,
    silent::Bool=false,
    )

    @assert length(eigVals) == length(eigVecs)

    @assert issorted(freqValues)

    broadeningFunc(x, standDev) = standDev ./ (x .^ 2 .+ standDev .^ 2)

    energyGs = minimum(eigVals)
    specFunc = 0 .* freqValues

    degenerateManifold = eigVals .≤ energyGs + degenTol
    if !silent
        println("Degeneracy = ", length(eigVals[degenerateManifold]), "; Range=[$(eigVals[degenerateManifold][1]), $(eigVals[degenerateManifold][end])]")
    end

    @showprogress for groundState in eigVecs[degenerateManifold]
        excitationCreate = probes["create"] * groundState
        excitationDestroy = probes["destroy"] * groundState
        excitationCreateBra = groundState' * probes["create"]
        excitationDestroyBra = groundState' * probes["destroy"]
        for index in eachindex(eigVals)
            excitedState = eigVecs[index]
            spectralWeights = [(excitationDestroyBra * excitedState) * (excitedState' * excitationCreate),
                               (excitedState' * excitationDestroy) * (excitationCreateBra * excitedState)
                              ]
            specFunc .+= spectralWeights[1] * broadeningFunc(freqValues .+ energyGs .- eigVals[index], standDev)
            specFunc .+= spectralWeights[2] * broadeningFunc(freqValues .- energyGs .+ eigVals[index], standDev)
        end
    end
    areaSpecFunc = sum(specFunc .* (maximum(freqValues) - minimum(freqValues[1])) / (length(freqValues) - 1))
    if areaSpecFunc > 1e-10 && normalise
        specFunc = specFunc ./ areaSpecFunc
    end
    return specFunc
end
export SpecFunc


"""
    SpecFunc(eigVals, eigVecs, probe, probeDiag, freqValues, standDev)

Extends SpecFunc above by allowing passing the spectrum in the form native to fermions.jl (dictionaries).

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
    standDev::Union{Vector{Float64}, Float64};
    degenTol::Float64=0.,
    normalise::Bool=true,
    silent::Bool=false,
)
    eigenStates = Vector{Float64}[] 
    for vector in eigVecs
        push!(eigenStates, ExpandIntoBasis(vector, basisStates))
    end
    probeMatrices = Dict{String,Matrix{Float64}}(name => zeros(length(basisStates), length(basisStates)) for name in keys(probes))
    for (name,probe) in collect(probes)
        probeMatrices[name] = OperatorMatrix(basisStates, probe)
    end

    return SpecFunc(eigVals, eigenStates, probeMatrices, freqValues, standDev;
                    normalise=normalise, degenTol=degenTol, silent=silent)
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
    standDev::Union{Vector{Float64}, Float64},
    symmetries::Vector{Char};
    degenTol::Float64=0.,
    normalise::Bool=true,
    silent::Bool=false,
)
    classifiedSpectrum, classifiedEnergies = ClassifyBasis(eigVecs, symmetries; energies=eigVals)
    groundStateEnergy = minimum(eigVals)
    minimalEigVecs = Dict{BitVector,Float64}[]
    minimalEigVals = Float64[]
    groundState = eigVecs[sortperm(eigVals)[1]]
    for groundState in eigVecs[eigVals .≤ groundStateEnergy  + abs(groundStateEnergy) * degenTol]
        push!(minimalEigVals, groundStateEnergy)
        push!(minimalEigVecs, groundState)

        # calculate sector of excited states
        for operator in ["create", "destroy"]
            excitedState = ApplyOperator(probes[operator], groundState)
            if !isempty(excitedState)
                sector = GetSector(excitedState, symmetries)
                append!(minimalEigVecs, classifiedSpectrum[sector])
                append!(minimalEigVals, classifiedEnergies[sector])
            end
        end
    end

    @assert groundStateEnergy == minimum(minimalEigVals)

    return SpecFunc(minimalEigVals, minimalEigVecs, probes, freqValues, basisStates, standDev;
                    normalise=normalise, degenTol=degenTol, silent=silent)
end
export SpecFunc


"""
    SpecFunc(eigVals, eigVecs, probe, probeDiag, freqValues, 
             standDev, symmetries, groundStateSector)

Extends SpecFunc() by allowing specific symmetry sectors from which to choose the ground state.

# Examples
```jldoctest
julia> SpecFunc(eigenVals, eigenStates, probe, probeDag, freqValues, 1e-2, 
                ['N']; groundStateSector=(3, 0))
```
"""
function SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probes::Dict{String,Vector{Tuple{String,Vector{Int64},Float64}}},
    freqValues::Vector{Float64},
    basisStates::Vector{Dict{BitVector,Float64}},
    standDev::Union{Vector{Float64}, Float64},
    symmetries::Vector{Char},
    groundStateSector::Union{Tuple{Int64}, Tuple{Int64,Int64}};
    degenTol::Float64=0.,
    normalise::Bool=true,
    silent::Bool=false,
)
    classifiedSpectrum, classifiedEnergies = ClassifyBasis(eigVecs, symmetries; energies=eigVals)
    minimalEigVecs = Dict{BitVector,Float64}[]
    minimalEigVals = Float64[]

    groundStateEnergy = minimum(classifiedEnergies[groundStateSector])

    groundStateCandidates = classifiedSpectrum[groundStateSector]
    for groundState in groundStateCandidates[classifiedEnergies[groundStateSector] .≤ groundStateEnergy  + abs(groundStateEnergy) * degenTol]

        # calculate sector of excited states
        push!(minimalEigVals, groundStateEnergy)
        push!(minimalEigVecs, groundState)
        for operator in ["create", "destroy"]
            excitedState = ApplyOperator(probes[operator], groundState)
            if !isempty(excitedState)
                sector = GetSector(excitedState, symmetries)
                append!(minimalEigVecs, classifiedSpectrum[sector])
                append!(minimalEigVals, classifiedEnergies[sector])
            end
        end
    end

    @assert groundStateEnergy == minimum(minimalEigVals)

    return SpecFunc(minimalEigVals, minimalEigVecs, probes, freqValues, basisStates, standDev;
                    normalise=normalise, degenTol=degenTol, silent=silent)
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
    standDev::Union{Vector{Float64}, Float64};
    degenTol::Float64=0.,
    normalise::Bool=true,
    silent::Bool=false,
)

    eigVecs = [collect(vec) for vec in eachcol(eigVecMatrix)]
    return SpecFunc(eigVals, eigVecs, probes, freqValues, standDev; 
                    normalise=normalise, degenTol=degenTol, silent=silent)
end
export SpecFunc
