function expandBasis(basisStates::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}}, 
        numAdditionalSites::Int64,
        energyVals::Dict{Tuple{Int64, Int64}, Vector{Float64}},
        retainSize::Int64;
        totOccReq::Vector{Int64}=[],
        totSzReq::Vector{Int64}=[]
    )
    @assert numAdditionalSites > 0
    newBasisStates = typeof(basisStates)()
    diagElements = Dict{Tuple{Int64, Int64}, Vector{Float64}}()
    for newBitCombination in collect(Iterators.product(repeat([(0,1),], 2 * numAdditionalSites)...))
        extraOcc = sum(newBitCombination)
        extraSz = sum(newBitCombination[1:2:end]) - sum(newBitCombination[2:2:end])

        Threads.@threads for ((totOcc, totSz), basisArr) in collect(basisStates)
            if (!isempty(totOccReq) && (totOcc + extraOcc) ∉ totOccReq) || (!isempty(totSzReq) && (totSz + extraSz) ∉ totSzReq)
                continue
            end
            newKey = (totOcc + extraOcc, totSz + extraSz)
            for (energy, basisStateDict) in zip(energyVals[(totOcc, totSz)], basisArr)
                expandedKeys = [vcat(key, collect(newBitCombination)) for key in keys(basisStateDict)]
                if newKey ∉ keys(newBasisStates)
                    newBasisStates[newKey] = []
                    diagElements[newKey] = []
                end
                push!(newBasisStates[newKey], Dict(zip(expandedKeys, values(basisStateDict))))
                push!(diagElements[newKey], energy)
            end
        end
    end
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


function iterativeDiagonaliser(
        hamiltonianFamily::Vector{Dict{Tuple{String,Vector{Int64}},Float64}},
        initBasis::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}},
        numStatesFamily::Vector{Int64},
        retainSize::Int64;
        tolerance::Float64
    )
    @assert length(hamiltonianFamily) == length(numStatesFamily)
    spectrumFamily = []
    basisStates = initBasis
    diagElements = Dict{Tuple{Int64, Int64}, Vector{Float64}}(k => zeros(length(v)) for (k, v) in basisStates)
    @showprogress for (i, hamiltonian) in enumerate(hamiltonianFamily)
        hamiltonianMatrix = operatorMatrix(basisStates, hamiltonian)
        for k in keys(hamiltonianMatrix)
            hamiltonianMatrix[k] .+= diagm(diagElements[k])
        end
        spectrum = getSpectrum(hamiltonianMatrix; maxNum=retainSize, tolerance=tolerance)
        push!(spectrumFamily, spectrum)
        @assert keys(basisStates) == keys(spectrum[1]) == keys(spectrum[2])
        basisStates = transformBasis(basisStates, spectrum[2])
        if i < length(numStatesFamily)
            basisStates, diagElements = expandBasis(basisStates, numStatesFamily[i+1] - numStatesFamily[i], spectrum[1], retainSize; totOccReq=[1+numStatesFamily[i+1]], totSzReq=[1, 0, -1])
        end
    end
    return spectrumFamily
end
