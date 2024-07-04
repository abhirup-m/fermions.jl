"""Expands the basis to accomodate new 1-particle states by tacking 
product states on to existing states. For eg. |10> + |01> -> (|10> + |01>)⊗|1>.
"""
function expandBasis(
        basisStates::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}}, 
        numAdditionalSites::Int64,
        energyVals::Dict{Tuple{Int64, Int64}, Vector{Float64}};
        occCriteria::Function = (x,y) -> true,
        magzCriteria::Function = x -> true
    )
    @assert numAdditionalSites > 0

    # for each new basis state, the variable 'diagElements' stores the energy 
    # of the old eigenstate from which the new basis state was derived. Serves
    # two purposes:(i) sorts the full set of new basis states so that we can 
    # take the M least energetic states (truncation), and (ii) these energy
    # values form the diagonal elements of the new Nth Hamiltonian.
    diagElements = Dict{Tuple{Int64, Int64}, Vector{Float64}}()

    numLevelsOld = length.(keys(collect(values(basisStates))[1][1]))[1]

    # each new basis states is associated with an old basis state from which 
    # it was derived. These old basis states are actually eigenstates of the
    # N-1-th hamiltonian in an iterative diagonalisation scheme, and are
    # themselves associated with energy eigenvalues.

    allowedSectors = Dict{Tuple{Int64, Int64}, Vector{Tuple{Tuple{Int64, Int64}, Vector{Int64}}}}()
    for (occ, magz) in keys(basisStates) 
        for newBitCombination in Iterators.product(repeat([(0,1),], 2 * numAdditionalSites)...)
            totalOcc = occ + sum(newBitCombination)
            totalMagz = magz + sum(newBitCombination[1:2:end]) - sum(newBitCombination[2:2:end])
            numLevelsNew = numLevelsOld + length(newBitCombination)
            @assert numLevelsNew % 2 == 0
            if !occCriteria(totalOcc, trunc(Int, numLevelsNew / 2)) || !magzCriteria(totalMagz)
                continue
            end
            if (totalOcc, totalMagz) ∉ keys(allowedSectors)
                allowedSectors[(totalOcc, totalMagz)] = []
            end
            push!(allowedSectors[(totalOcc, totalMagz)], ((occ, magz), collect(newBitCombination)))
            if (totalOcc, totalMagz) ∉ keys(diagElements)
                diagElements[(totalOcc, totalMagz)] = []
            end
            diagElements[(totalOcc, totalMagz)] = [diagElements[(totalOcc, totalMagz)]; energyVals[(occ, magz)]]
        end
    end

    # stores the new expanded basis
    newBasisStates = Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}}(k => [] for k in keys(allowedSectors))

    @showprogress desc="expand basis" Threads.@threads for (newSector, possibilities) in collect(allowedSectors)

        newBasisStates[newSector] = fetch.([Threads.@spawn ((x,y) -> Dict(vcat(k, x) => v for (k,v) in y))(newBitCombination, basisStateDict)
                                            for (oldSector, newBitCombination) in possibilities 
                                            for (energy, basisStateDict) in zip(energyVals[oldSector],basisStates[oldSector])
                                           ])

        """
            # loop over the computational basis components of the old state
            for (energy, basisStateDict) in zip(energyVals[oldSector], basisStates[oldSector])
                # form new basis states by joining the new and the old 
                # computational states: |10> -><- |0> = |100>
                @time expandedBasisStates = [vcat(key, newBitCombination)
                                       for key in keys(basisStateDict)]

                @time push!(newBasisStates[newSector], Dict(zip(expandedBasisStates, values(basisStateDict))))
                @time push!(diagElements[newSector], energy)
            end
        end
        """
        if isempty(newBasisStates[newSector])
            delete!(newBasisStates, newSector)
            delete!(diagElements, newSector)
        end
    end

    return newBasisStates, diagElements
end


"""Main function for iterative diagonalisation. Gives the approximate low-energy
spectrum of a hamiltonian through the following algorithm: first diagonalises a
Hamiltonian with few degrees of freedom, retains a fixed Ns number of low-energy
eigenstates, writes the Hamiltonian for a larger number of degrees of freedom in
the basis of these Ns states, then diagonalises it, etc.
"""
function iterativeDiagonaliser(
        hamiltonianFamily::Vector{Dict{Tuple{String,Vector{Int64}},Float64}},
        initBasis::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}},
        numStatesFamily::Vector{Int64},
        retainSize::Int64;
        occCriteria::Function = (x,y) -> true,
        magzCriteria::Function = x -> true,
        keepfrom::String="IR",
        tolerance::Float64=1e-12
    )
    @assert length(hamiltonianFamily) == length(numStatesFamily)

    # stores flow of the spectrum
    spectrumFamily = []

    # stores the flowing basis states
    basisStates = initBasis

    # stores the current eigenvalues
    diagElements = Dict{Tuple{Int64, Int64}, Vector{Float64}}(k => zeros(length(v)) for (k, v) in basisStates)

    # loop over the Hamiltonians H_N.
    @showprogress desc="iterative diag." for (i, hamiltonian) in enumerate(hamiltonianFamily)

        # at Nth step on bringing sites 'q' into the Hamiltonian,
        # full hamiltonian H_N is H_N = H_{N-1} + H_{N-1,q} + H_q,
        # where H_q is the terms acting purely on q, while H_{N-1,q}
        # is the set of terms that couple 'q' with the N-1 states.

        # matrix for the H_{N-1,q} + H_q part.
        hamiltonianMatrix = operatorMatrix(basisStates, hamiltonian)
        # Adds the matrix elements for H_{N-1}. Since the present basis
        # is the eigenbasis of H_{N-1}, these matrix elements are the diagonal ones.
        for k in keys(hamiltonianMatrix)
            hamiltonianMatrix[k] .+= diagm(diagElements[k])
        end

        # diagonalise this matrix, keeping only 'retainSize' 
        # number of the new eigenstates.
        spectrum = getSpectrum(basisStates, hamiltonianMatrix; maxNum=retainSize, keepfrom=keepfrom, tolerance=tolerance)

        push!(spectrumFamily, spectrum)
        @assert keys(basisStates) == keys(spectrum[1]) == keys(spectrum[2])
        # create basis out of the eigenstates
        basisStates = spectrum[2]

        if i < length(numStatesFamily)
            # expand the new basis to accomodate the new sites for the next step.
            basisStates, diagElements = expandBasis(basisStates, numStatesFamily[i+1] - numStatesFamily[i], spectrum[1]; 
                                                          occCriteria=occCriteria, magzCriteria=magzCriteria)
        end
    end
    return spectrumFamily
end
