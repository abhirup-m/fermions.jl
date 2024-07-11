using IterTools
include("base.jl")
include("eigen.jl")

"""Expands the basis to accomodate new 1-particle states by tacking 
product states on to existing states. For eg. |10> + |01> -> (|10> + |01>)⊗|1>.
"""
function ExpandBasis(
    basis::Vector{Dict{BitVector,Float64}},
    sector::Tuple{Int64,Int64},
    eigvals::Vector{Float64},
    numNew::Int64;
    newClassifiedBasis::Dict{Tuple{Int64,Int64},Vector{Dict{BitVector,Float64}}}=Dict{Tuple{Int64,Int64},Vector{Dict{BitVector,Float64}}}(),
    newDiagElementsClassified::Dict{Tuple{Int64,Int64},Vector{Float64}}=Dict{Tuple{Int64,Int64},Vector{Float64}}(),
)
    @assert length(basis) == length(eigvals)
    newPairs = []
    for newComb in collect(IterTools.product(repeat([[0, 1]], 2 * numNew)...))
        newSector = sector .+ (sum(newComb),
            sum(newComb[1:2:end]) - sum(newComb[2:2:end]))
        if newSector ∉ keys(newClassifiedBasis)
            newClassifiedBasis[newSector] = []
            newDiagElementsClassified[newSector] = []
        end
        push!(newPairs, (newComb, newSector))
    end
    Threads.@threads for (newComb, newSector) in newPairs
        append!(newClassifiedBasis[newSector], [Dict([vcat(key, collect(newComb)) for key in keys(stateDict)] .=> values(stateDict)) for stateDict in basis])
        append!(newDiagElementsClassified[newSector], eigvals)
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
    occupSectors::Function=(x, N) -> true
)
    @assert length(numSitesFlow) == length(hamFlow)

    classifiedBasis = ClassifyBasis(basisStates)
    diagElementsClassified = Dict(k => zeros(length(v)) for (k, v) in classifiedBasis)

    numDiffFlow = [n2 - n1 for (n1, n2) in zip(numSitesFlow[1:end-1], numSitesFlow[2:end])]
    eigvalFlow = [Dict() for _ in hamFlow]
    eigvecFlow = [Dict() for _ in hamFlow]

    for (i, hamiltonian) in enumerate(hamFlow)
        newClassifiedBasis = typeof(classifiedBasis)()
        newDiagElementsClassified = Dict{keytype(classifiedBasis),Vector{Float64}}()
        spectrum = fetch.([Threads.@spawn Spectrum(hamiltonian, basis, diagElements=diagElementsClassified[sector]) for (sector, basis) in classifiedBasis])
        for (sector, spectrumSector) in zip(keys(classifiedBasis), spectrum)
            eigvalFlow[i][sector], eigvecFlow[i][sector] = spectrumSector
        end
        if i == length(hamFlow)
            continue
        end
        for sector in keys(classifiedBasis)
            newClassifiedBasis, newDiagElementsClassified = ExpandBasis(eigvecFlow[i][sector], sector,
                eigvalFlow[i][sector], numDiffFlow[i];
                newClassifiedBasis=newClassifiedBasis, newDiagElementsClassified=newDiagElementsClassified)
        end

        if i == length(hamFlow)
            continue
        end
        classifiedBasis, diagElementsClassified = TruncateBasis(newClassifiedBasis, newDiagElementsClassified, retainSize)
    end
    return eigvalFlow, eigvecFlow
end
