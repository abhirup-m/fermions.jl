using Serialization, Random, LinearAlgebra

"""Expands the basis to accomodate new 1-particle states by tacking 
product states on to existing states. For eg. |10> + |01> -> (|10> + |01>)⊗|1>.
"""
function ExpandBasis(
    basis::Vector{Dict{BitVector,Float64}},
    sector::Tuple{Int64,Int64},
    eigvals::Vector{Float64},
    numNew::Int64;
    totOccCriteria::Function=o -> true,
    magzCriteria::Function=m -> true,
    newClassifiedBasis::Dict{Tuple{Int64,Int64},Vector{Dict{BitVector,Float64}}}=Dict{Tuple{Int64,Int64},Vector{Dict{BitVector,Float64}}}(),
    newDiagElementsClassified::Dict{Tuple{Int64,Int64},Vector{Float64}}=Dict{Tuple{Int64,Int64},Vector{Float64}}(),
)
    @assert length(basis) == length(eigvals)
    newPairs = Tuple{NTuple{2 * numNew,Int64},Tuple{Int64,Int64}}[]
    for newComb in collect(Iterators.product(repeat([[0, 1]], 2 * numNew)...))
        newSector = sector .+ (sum(newComb),
            sum(newComb[1:2:end]) - sum(newComb[2:2:end]))
        if !(totOccCriteria(newSector[1]) && magzCriteria(newSector[2]))
            continue
        end
        if newSector ∉ keys(newClassifiedBasis)
            newClassifiedBasis[newSector] = []
            newDiagElementsClassified[newSector] = []
        end
        push!(newPairs, (newComb, newSector))
    end
    Threads.@threads for (newComb, newSector) in newPairs
        append!(newDiagElementsClassified[newSector], eigvals)
        for stateDict in basis
            newKeyValPairs = [(vcat(key, collect(newComb)), val) for (key, val) in stateDict]
            push!(newClassifiedBasis[newSector], Dict(newKeyValPairs))
        end
    end
    return newClassifiedBasis, newDiagElementsClassified
end


function TruncateBasis(newClassifiedBasis, newDiagElementsClassified, retainSize)
    classifiedBasis = typeof(newClassifiedBasis)()
    diagElementsClassified = Dict{keytype(classifiedBasis),Vector{Float64}}()
    for (sector, basis) in newClassifiedBasis
        sortperms = sortperm(newDiagElementsClassified[sector])
        classifiedBasis[sector] = basis[sortperms][1:minimum((retainSize, length(basis)))]
        diagElementsClassified[sector] = sort(newDiagElementsClassified[sector])[1:minimum((retainSize, length(newDiagElementsClassified[sector])))]
    end
    return classifiedBasis, diagElementsClassified
end


"""Main function for iterative diagonalisation. Gives the approximate low-energy
spectrum of a hamiltonian through the following algorithm: first diagonalises a
Hamiltonian with few degrees of freedom, retains a fixed Ns number of low-energy
eigenstates, writes the Hamiltonian for a larger number of degrees of freedom in
the basis of these Ns states, then diagonalises it, etc.
"""
function IterDiag(
    hamFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    initBasis::Vector{Dict{BitVector,Float64}},
    numSitesFlow::Vector{Int64},
    maxSize::Int64;
    showProgress::Bool=false,
    totOccCriteria::Function=(o, N) -> true,
    magzCriteria::Function=m -> true,
)
    @assert length(numSitesFlow) == length(hamFlow)
    fundMatrices = Dict{Char, Matrix{Float64}}(type => OperatorMatrix(BasisStates(1), [(string(type), [1], 1.0)]) for type in ('+', '-', 'n', 'h'))
    sigmaz = diagm([1, -1])
    basicMatrices = Dict{Tuple{Char, Int64}, Matrix{Float64}}((type, site) => OperatorMatrix(basisStates, [(string(type), [site], 1.0)]) for site in 1:numSitesFlow[1] for type in ('+', '-', 'n', 'h'))
    transformMatrices = []
    eigvalData = []

    hamMatrix = TensorProduct(hamFlow[1], basicMatrices)
    F = eigen(hamMatrix)
    push!(transformMatrices, F.vectors[:, 1:maxSize])
    push!(eigvalData, F.values[1:maxSize])
    numSites = numSitesFlow[1]

    for (i, addedSites) in enumerate(numSitesFlow[2:end] .- numSitesFlow[1:end-1])
        for (k, v) in basicMatrices
            basicMatrices[k] = kron(transformMatrices[end]' * v * transformMatrices[end], I(2))
        end
        
        padSize = Int(transformMatrices[end] |> size |> minimum |> log2 |> ceil)
        identityRest = ifelse(padSize > 1, kron(fill(sigmaz, padSize)...), sigmaz)
        identityRestTransf = transformMatrices[end]' * identityRest * transformMatrices[end]
        display(identityRestTransf)

        for type in ('+', '-', 'n', 'h')
            display(fundMatrices[type])
            basicMatrices[(type, numSites+1)] = kron(identityRestTransf, fundMatrices[type])
        end

        display(basicMatrices[('+', 3)])
        display(basicMatrices[('-', 3)])
        hamMatrix = diagm(repeat(eigvalData[end], inner=2 * addedSites)) + TensorProduct(hamFlow[i+1], basicMatrices)
        display(hamMatrix)
        F = eigen(hamMatrix)
        push!(transformMatrices, F.vectors[:, 1:maxSize])
        push!(eigvalData, F.values)
        numSites += addedSites
    end
    return eigvalData, transformMatrices
end


function IterSpecFunc(groundSectors, spectrumData, probe, probeDag, freqArray, broadening, symmetries)
    specfuncFlow = [0 .* freqArray for _ in spectrumData]

    for (i, (sector, (eigVals, eigVecs))) in enumerate(zip(groundSectors, spectrumData))
        specfuncFlow[i] .+= SpecFunc(eigVals, eigVecs, probe, probeDag, freqArray, broadening, sector, symmetries)
    end
    specfuncFlow[end] ./= sum(specfuncFlow[end]) * abs(freqArray[2] - freqArray[1])
    return specfuncFlow, freqArray
end


function TensorProduct(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        basicMatrices::Dict{Tuple{Char, Int64}, Matrix{Float64}},
    )
    totalMat = sum([opStrength * prod([basicMatrices[pair] for pair in zip(opType, opMembers)]) for (opType, opMembers, opStrength) in operator])
end
