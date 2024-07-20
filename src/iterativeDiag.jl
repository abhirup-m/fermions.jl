"""Expands the basis to accomodate new 1-particle states by tacking 
product states on to existing states. For eg. |10> + |01> -> (|10> + |01>)⊗|1>.
"""
function ExpandBasis(
    basis::Vector{Dict{BitVector,Float64}},
    sector::Tuple{Int64,Int64},
    eigvals::Vector{Float64},
    retainSize::Int64,
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
            minVal = abs(sort(collect(values(stateDict)), rev=true)[minimum((retainSize, length(stateDict)))])
            newKeyValPairs = [(vcat(key, collect(newComb)), val) for (key, val) in stateDict if abs(val) >= minVal]
            push!(newClassifiedBasis[newSector], Dict(newKeyValPairs))
        end
    end
    return newClassifiedBasis, newDiagElementsClassified
end


function ClassifyBasis(basisStates::Vector{Dict{BitVector,Float64}})
    classifiedBasis = Dict{Tuple{Int64,Int64},typeof(basisStates)}()

    # Multhreading doesn't help here, already tried!
    for stateDict in basisStates
        bstate = collect(keys(stateDict))[1]
        totOcc = sum(bstate)
        totMagzArr = sum(bstate[1:2:end]) - sum(bstate[2:2:end])
        if (totOcc, totMagzArr) ∉ keys(classifiedBasis)
            classifiedBasis[(totOcc, totMagzArr)] = []
        end
        push!(classifiedBasis[(totOcc, totMagzArr)], stateDict)
    end

    return classifiedBasis
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
    basisStates::Vector{Dict{BitVector,Float64}},
    numSitesFlow::Vector{Int64},
    retainSize::Int64;
    totOccCriteria::Function=(o, N) -> true,
    magzCriteria::Function=m -> true,
)
    @assert length(numSitesFlow) == length(hamFlow)

    classifiedBasis = ClassifyBasis(basisStates)
    diagElementsClassified = Dict(k => zeros(length(v)) for (k, v) in classifiedBasis)

    numDiffFlow = [n2 - n1 for (n1, n2) in zip(numSitesFlow[1:end-1], numSitesFlow[2:end])]
    spectrumFlow = []

    @showprogress for (i, hamiltonian) in enumerate(hamFlow)
        push!(spectrumFlow, [])
        newClassifiedBasis = typeof(classifiedBasis)()
        newDiagElementsClassified = Dict{keytype(classifiedBasis),Vector{Float64}}()
        spectrumFlow[end] = Dict(keys(classifiedBasis) .=> fetch.([Threads.@spawn Spectrum(hamiltonian, basis, diagElements=diagElementsClassified[sector])
                                                                   for (sector, basis) in classifiedBasis]))
        if i == length(hamFlow)
            continue
        end
        for sector in keys(classifiedBasis)
            newClassifiedBasis, newDiagElementsClassified = ExpandBasis(spectrumFlow[end][sector][2], sector,
                spectrumFlow[end][sector][1], retainSize, numDiffFlow[i]; totOccCriteria=o -> totOccCriteria(o, numSitesFlow[i+1]),
                magzCriteria=magzCriteria, newClassifiedBasis=newClassifiedBasis,
                newDiagElementsClassified=newDiagElementsClassified)
        end

        if i == length(hamFlow)
            continue
        end
        classifiedBasis, diagElementsClassified = TruncateBasis(newClassifiedBasis, newDiagElementsClassified, retainSize)
    end
    return spectrumFlow
end
