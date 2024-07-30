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
    basisStates::Vector{Dict{BitVector,Float64}},
    numSitesFlow::Vector{Int64},
    retainSizeY::Int64;
    showProgress::Bool=false,
    totOccCriteria::Function=(o, N) -> true,
    magzCriteria::Function=m -> true,
)
    @assert length(numSitesFlow) == length(hamFlow)

    classifiedBasis = ClassifyBasis(basisStates, ['N', 'Z'])
    diagElementsClassified = Dict(k => zeros(length(v)) for (k, v) in classifiedBasis)

    numDiffFlow = [n2 - n1 for (n1, n2) in zip(numSitesFlow[1:end-1], numSitesFlow[2:end])]
    eigValData = Dict{Tuple{Int64, Int64}, Vector{Float64}}[]
    eigVecData = Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}}[]

    @showprogress enabled=showProgress for (i, hamiltonian) in enumerate(hamFlow)
        push!(eigValData, Dict())
        push!(eigVecData, Dict())
        eigValData[end], eigVecData[end] = Spectrum(hamiltonian, classifiedBasis, diagElements=diagElementsClassified)
        if i == length(hamFlow)
            continue
        end

        newClassifiedBasis = typeof(classifiedBasis)()
        newDiagElementsClassified = Dict{keytype(classifiedBasis),Vector{Float64}}()
        for sector in keys(classifiedBasis)
            newClassifiedBasis, newDiagElementsClassified = ExpandBasis(eigVecData[end][sector], sector,
                eigValData[end][sector], numDiffFlow[i]; totOccCriteria=o -> totOccCriteria(o, numSitesFlow[i+1]),
                magzCriteria=magzCriteria, newClassifiedBasis=newClassifiedBasis,
                newDiagElementsClassified=newDiagElementsClassified)
        end

        if i == length(hamFlow)
            continue
        end
        classifiedBasis, diagElementsClassified = TruncateBasis(newClassifiedBasis, newDiagElementsClassified, retainSizeY)
    end
    dirStructure = joinpath("data")
    mkpath(dirStructure)
    savePathE = joinpath(dirStructure, string(now()) |> x->replace(x,":"=>"-","."=>"-")) * "-E"
    savePathX = joinpath(dirStructure, string(now()) |> x->replace(x,":"=>"-","."=>"-")) * "-X"
    serialize(savePathE, eigValData)
    serialize(savePathX, eigVecData)
    return eigValData, eigVecData
end


function IterSpecFunc(groundSectors, spectrumData, probe, probeDag, freqArray, broadening, symmetries)
    specfuncFlow = [0 .* freqArray for _ in spectrumData]

    for (i, (sector, (eigVals, eigVecs))) in enumerate(zip(groundSectors, spectrumData))
        specfuncFlow[i] .+= SpecFunc(eigVals, eigVecs, probe, probeDag, freqArray, broadening, sector, symmetries)
    end
    specfuncFlow[end] ./= sum(specfuncFlow[end]) * abs(freqArray[2] - freqArray[1])
    return specfuncFlow, freqArray
end
