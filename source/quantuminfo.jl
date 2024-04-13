using LinearAlgebra

function ReducedDensityMatrix(genState::Dict{String, Float64}, partiesRemain)
    if length(partiesRemain) == length(genState)
        redDenMat = values(genState) * collect(values(genState)')
        return redDenMat
    end
    getsubstring(s, indices) = join([ch for (i, ch) in enumerate(s) if i in indices])

    remainBasis = BasisStates(length(partiesRemain))

    partiesTraced = setdiff(1:length(collect(keys(genState))[1]), partiesRemain)
    stateTracedBasis = Dict()
    for (state, coeff) in genState
        labelRemain = getsubstring(state, partiesRemain)
        labelTraced = getsubstring(state, partiesTraced)
        if ! haskey(stateTracedBasis, labelTraced) 
            stateTracedBasis[labelTraced] = Dict([(state, 0.) for state in remainBasis])
    	end
        stateTracedBasis[labelTraced][labelRemain] += coeff
    end

    redDenMat = (+)([collect(values(stateTracedBasis[labelTraced])) * collect(values(stateTracedBasis[labelTraced]))'
		     for labelTraced in keys(stateTracedBasis)]...)

    return redDenMat
end

function EntanglementEntropy(genState, parties)
    redDenMat = ReducedDensityMatrix(genState, parties)
    eigenvalues = eigvals(Hermitian(redDenMat))
    nonzero_eigvals = eigenvalues[eigenvalues .> 0]
    entEntropy = -sum(nonzero_eigvals .* log.(nonzero_eigvals))
    return entEntropy
end

function MutualInfo(genState, parties)

    S_A = EntanglementEntropy(genState, parties[1])
    S_B = EntanglementEntropy(genState, parties[2])
    S_AB = EntanglementEntropy(genState, [parties[1]; parties[2]])

    return S_A + S_B - S_AB
end


function mapToFunction(measure)
    if measure == "VNE"
        return EntanglementEntropy
    elseif measure == "I2"
        return MutualInfo
    else
        return
    end
end