"""
Returns a set of basis states in the form
[[0,0,...,0],
 [0,0,...,1],
       .
       .
       .
 [1,1,...,1]
],
where each inner vector represents a configuration of the electronic fock states.
Optionally accepts a parameter totOccupancy that restricts the basis states to only
those with the specified total occupancy, and a parameter totSpin that does the same
but for a specified total magnetisation.
"""


function roundTo(val, tolerance)
    return round(val, digits=trunc(Int, -log10(tolerance)))
end


function multiplyDict(dicts, multipliers)
    @assert length(dicts) == length(multipliers)
    return mergewith(+, fetch.([Threads.@spawn ((dict, m) -> Dict(k => v * m for (k, v) in dict))(dict, m) for (dict, m) in zip(dicts, multipliers) if m != 0])...)
end


function BasisStates(numLevels::Int64; magzCriteria::Function=x -> true, occCriteria::Function=x -> true, localConstraint::Function=x -> true)
    @assert numLevels > 0 && numLevels % 2 == 0

    # create dictionary to store configurations classified by their total occupancies.
    # For eg, basisStates = {0: [[0, 0]], 1: [[1, 0], [0, 1]], 2: [[1, 1]]}}
    basisStates = Dict{Tuple{Int64,Int64},Vector{Dict{BitVector, Float64}}}()

    # generate all possible configs, by running integers from 0 to 2^N-1, and converting
    # them to binary forms
    allConfigs = digits.(0:2^numLevels-1, base=2, pad=numLevels) |> reverse

    for config in allConfigs

        # check which states have the correct occupancy and total Sz
        occupancy = sum(config)
        magz = sum(config[1:2:end]) - sum(config[2:2:end])
        if !occCriteria(occupancy) || !magzCriteria(magz) || !localConstraint(config)
            continue
        end

        if !((occupancy, magz) in keys(basisStates))
            basisStates[(occupancy, magz)] = []
        end
        if length(config) > 0
            push!(basisStates[(occupancy, magz)], Dict(config => 1.0))
        end
    end

    return basisStates
end


function TransformBit(qubit::Bool, operator::Char)
    @assert operator in ('n', 'h', '+', '-')
    if operator == 'n'
        return qubit, qubit
    elseif operator == 'h'
        return qubit, 1 - qubit
    elseif (operator == '+' && qubit == 0) || (operator == '-' && qubit == 1)
        return 1 - qubit, 1
    else
        return qubit, 0
    end
end


function applyOperatorOnState(stateDict::Dict{BitVector,Float64}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    @assert maximum([maximum(opMembers) for (_, opMembers) in keys(operatorList)]) ≤ length(collect(keys(stateDict))[1])

    # define a dictionary for the final state obtained after applying 
    # the operators on all basis states
    completeOutputState = Dict{BitVector,Float64}()
    # loop over all operator tuples within operatorList
    for ((opType, opMembers), opStrength) in pairs(operatorList)

        for (state, coefficient) in stateDict

            newCoefficient = coefficient
            newState = copy(state)

            # for each basis state, obtain a modified state after applying the operator tuple
            for (siteIndex, operator) in zip(reverse(opMembers), reverse(opType))
                newQubit, factor = TransformBit(newState[siteIndex], operator)
                if factor == 0
                    newCoefficient = 0
                    break
                end
                # calculate the fermionic exchange sign by counting the number of
                # occupied states the operator has to "hop" over
                exchangeSign = operator in ['+', '-'] ? (-1)^sum(newState[1:siteIndex-1]) : 1
                @inbounds newState[siteIndex] = newQubit
                newCoefficient *= exchangeSign * factor
            end
            if newCoefficient != 0
                if newState in keys(completeOutputState)
                    @inbounds completeOutputState[newState] += opStrength * newCoefficient
                else
                    @inbounds completeOutputState[newState] = opStrength * newCoefficient
                end
            end
        end
    end
    return completeOutputState
end


"""Creates matrix representation for the provided abstract operator definition.
"""
function operatorMatrix(basisStates::Dict{Tuple{Int64,Int64}, Vector{Dict{BitVector, Float64}}}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64}; tolerance=1e-14)
    operatorFullMatrix = Dict(key => zeros(length(value), length(value)) for (key, value) in basisStates)

    # loop over symmetry sectors of the basis
    Threads.@threads for (key, bstates) in collect(basisStates)

        # loop over basis states within the symmetry sector.
        # These states will act as incoming states in terms
        # of the operator (hence the column index 'col').
        Threads.@threads for (col, colState) in collect(enumerate(bstates))

            # get action of operator on incoming basis state |ψ>
            newStateDict = applyOperatorOnState(colState, operatorList)

            # again loop over the basis states, this term for the outgoing states <ϕ|
            for (row, rowState) in enumerate(bstates)
                if operatorFullMatrix[key][row, col] != 0
                    continue
                end

                # calculate overlap between the outgoing state <ϕ| and the 
                # modified incoming state O|ψ>.
                overlap = 0
                for key in keys(newStateDict)
                    if key in keys(rowState)
                    overlap += newStateDict[key] * rowState[key]
                    end
                end
                operatorFullMatrix[key][row, col] = abs(overlap) > tolerance ? overlap : 0
                operatorFullMatrix[key][col, row] = operatorFullMatrix[key][row, col]'
            end
        end
    end
    return operatorFullMatrix
end


function prettyPrint(state::BitVector)
    output = ""
    @assert length(state) % 2 == 0
    for i in 1:2:length(state)
        if state[[i, i + 1]] == [0, 0]
            output = output * "0"
        elseif state[[i, i + 1]] == [1, 0]
            output = output * "↑"
        elseif state[[i, i + 1]] == [0, 1]
            output = output * "↓"
        else
            output = output * "2"
        end
    end
    return output
end


function getSpectrum(
        basisStates::Dict{Tuple{Int64,Int64},Vector{Dict{BitVector, Float64}}},
        hamiltonian::Dict{Tuple{Int64,Int64},Matrix{Float64}}; 
        tolerance::Float64=1e-14, maxNum::Int64=0, keepfrom="IR",
        progressEnabled::Bool=false
    )
    for (k, v) in hamiltonian
        if !ishermitian(v) 
            if abs(maximum(v .- v')) < tolerance
                hamiltonian[k] = 0.5 .* (v .+ v')
            else
                @assert false "Not hermitian."
            end
        end
    end

    eigvals = Dict{Tuple{Int64,Int64},Vector{Float64}}(index => [] for index in keys(hamiltonian))
    eigvecs = Dict{Tuple{Int64,Int64},Vector{Dict{BitVector, Float64}}}(index => [] for index in keys(hamiltonian))
    pbar = nothing
    if progressEnabled
        pbar = Progress(length(hamiltonian))
    end

    for (index, mat) in collect(hamiltonian)
        F = eigen(Hermitian(mat))
        maxNumSector = ifelse(maxNum == 0, length(F.values), minimum((length(F.values), maxNum)))
        sliceRange = keepfrom == "IR" ? (1:maxNumSector) : ((length(F.values) - maxNumSector + 1):length(F.values))
        @inbounds eigvals[index] = roundTo.(sort(F.values)[1:maxNumSector], tolerance)
        @inbounds eigvecs[index] = fetch.([Threads.@spawn (v -> multiplyDict(basisStates[index], v))(roundTo.(F.vectors[:, i], tolerance))
                                           for i in sortperm(F.values)[sliceRange]])
    end
    if !isnothing(pbar)
        finish!(pbar)
    end
    return eigvals, eigvecs
end


function StateOverlap(state1::Dict{BitVector, Float64}, state2::Dict{BitVector, Float64})
    shortState, longState = ifelse(length(state1) < length(state2), (state1, state2), (state2, state1))
    overlap = sum(fetch.([Threads.@spawn (k -> state1[k] * state2[k])(key) for key in keys(shortState) if haskey(longState, key)]))
    return overlap
end
