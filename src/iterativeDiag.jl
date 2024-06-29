"""Expands the basis to accomodate new 1-particle states by tacking 
product states on to existing states. For eg. |10> + |01> -> (|10> + |01>)⊗|1>.
"""
function expandBasis(basisStates::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}}, 
        numAdditionalSites::Int64,
        energyVals::Dict{Tuple{Int64, Int64}, Vector{Float64}},
        retainSize::Int64;
        totOccReq::Vector{Int64}=[],
        totSzReq::Vector{Int64}=[]
    )
    @assert numAdditionalSites > 0

    # stores the new expanded basis
    newBasisStates = typeof(basisStates)()

    # each new basis states is associated with an old basis state from which 
    # it was derived. These old basis states are actually eigenstates of the
    # N-1-th hamiltonian in an iterative diagonalisation scheme, and are
    # themselves associated with energy eigenvalues.

    # for each new basis state, the variable 'diagElements' stores the energy 
    # of the old eigenstate from which the new basis state was derived. Serves
    # two purposes:(i) sorts the full set of new basis states so that we can 
    # take the M least energetic states (truncation), and (ii) these energy
    # values form the diagonal elements of the new Nth Hamiltonian.
    diagElements = Dict{Tuple{Int64, Int64}, Vector{Float64}}()

    # loop over all combinations of the new states that are being added. if only
    # one states is being added, then 0 and 1. if two, then 00, 01, 10 and 11.
    for newBitCombination in collect(Iterators.product(repeat([(0,1),], 2 * numAdditionalSites)...))

        # check the occupancy and Sz associated with these combinations
        extraOcc = sum(newBitCombination)
        extraSz = sum(newBitCombination[1:2:end]) - sum(newBitCombination[2:2:end])

        # loop over the existing old basis states
        Threads.@threads for ((totOcc, totSz), basisArr) in collect(basisStates)

            # only consider further if the total new occupancy and total new Sz
            # (after attaching the extra states with the old state) fall in the 
            # range that are required (the required range is passed as function arg.)
            if (!isempty(totOccReq) && (totOcc + extraOcc) ∉ totOccReq) || (!isempty(totSzReq) && (totSz + extraSz) ∉ totSzReq)
                continue
            end

            newKey = (totOcc + extraOcc, totSz + extraSz)

            # loop over the computational basis components of the old state
            for (energy, basisStateDict) in zip(energyVals[(totOcc, totSz)], basisArr)

                # form new basis states by joining the new and the old 
                # computational states: |10> -><- |0> = |100>
                expandedBasisStates = [vcat(key, collect(newBitCombination)) for key in keys(basisStateDict)]

                if newKey ∉ keys(newBasisStates)
                    newBasisStates[newKey] = []
                    diagElements[newKey] = []
                end

                push!(newBasisStates[newKey], Dict(zip(expandedKeys, values(basisStateDict))))
                push!(diagElements[newKey], energy)
            end
        end
    end

    # retain only up to 'retainSize' number of states.
    # keep the ones with lowest energy in the old basis.
    Threads.@threads for (k, statesArr) in collect(newBasisStates)
        keepUpto = minimum((length(newBasisStates[k]), retainSize))
        newBasisStates[k] = newBasisStates[k][sortperm(diagElements[k])][1:keepUpto]
        diagElements[k] = sort(diagElements[k])[1:keepUpto]
    end
    return newBasisStates, diagElements
end


"""Performs a rotation of the given basis according to the given transformation. 
Can be used to, for eg., transfer into an eigen basis."""
function transformBasis(basisStates::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}},
        transformation::Dict{Tuple{Int64, Int64}, Vector{Vector{Float64}}})

    # check that the symmetry sectors of the transformation 
    # are part of the old basis
    @assert issubset(keys(transformation), keys(basisStates))

    # stores the new rotated basis
    newBasisStates = Dict{keytype(basisStates), valtype(basisStates)}(k => [Dict() for vi in v] for (k, v) in transformation)

    # loop over symmetry sectors of the transformation
    Threads.@threads for (sector, vectors) in collect(transformation)

        # loop over the eigen vectors of the transformation
        # (for eg. eigenvector = [-0.5, 0.5, 0.5, 0.5])
        Threads.@threads for (i, eigenvector) in collect(enumerate(vectors))

            # loop over the v_i^j
            for (j, multiplier) in enumerate(eigenvector)
                if multiplier == 0
                    continue
                end

                # create the new basis states: {|1>: v_i^1 * b_{1, old}, |2>: v_i^2 * b_{2, old} ... }
                multipliedState = Dict(q => multiplier * b for (q, b) in basisStates[sector][j])
                mergewith!(+, newBasisStates[sector][i], multipliedState)
            end
            @assert !isempty(newBasisStates[sector][i])
        end
    end
    return newBasisStates
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
        tolerance::Float64
    )
    @assert length(hamiltonianFamily) == length(numStatesFamily)

    # stores flow of the spectrum
    spectrumFamily = []

    # stores the flowing basis states
    basisStates = initBasis

    # stores the current eigenvalues
    diagElements = Dict{Tuple{Int64, Int64}, Vector{Float64}}(k => zeros(length(v)) for (k, v) in basisStates)

    # loop over the Hamiltonians H_N.
    @showprogress for (i, hamiltonian) in enumerate(hamiltonianFamily)

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
        spectrum = getSpectrum(hamiltonianMatrix; maxNum=retainSize, tolerance=tolerance)
        push!(spectrumFamily, spectrum)
        @assert keys(basisStates) == keys(spectrum[1]) == keys(spectrum[2])

        # create basis out of the eigenstates
        basisStates = transformBasis(basisStates, spectrum[2])
        if i < length(numStatesFamily)
            # expand the new basis to accomodate the new sites for the next step.
            basisStates, diagElements = expandBasis(basisStates, numStatesFamily[i+1] - numStatesFamily[i], spectrum[1], retainSize; totOccReq=[1+numStatesFamily[i+1]], totSzReq=[1, 0, -1])
        end
    end
    return spectrumFamily
end
