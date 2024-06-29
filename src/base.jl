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
function BasisStates(numLevels::Int64; totSz::Union{Int64,Vector{Int64},Nothing}=nothing, totOccupancy::Union{Int64,Vector{Int64},Nothing}=nothing, localOccupancy::Tuple{Vector{Int64}, Int64}=([1,2], -1))
    @assert numLevels > 0 && numLevels % 2 == 0

    # if totOccupancy is not set, all occupancies from 0 to N are allowed,
    # otherwise use only the provided totOccupancy
    if isnothing(totOccupancy)
        allowedOccupancies = 0:numLevels
    elseif typeof(totOccupancy) == Vector{Int64}
        @assert all(0 .<= totOccupancy .<= numLevels)
        allowedOccupancies = totOccupancy
    else
        @assert all(0 <= totOccupancy <= numLevels)
        allowedOccupancies = [totOccupancy]
    end

    # same for total Sz
    numSites = trunc(Int, numLevels / 2)
    if isnothing(totSz)
        allowedSz = -numSites:1:numSites
    elseif typeof(totSz) == Vector{Int64}
        @assert all(-numSites .<= totSz .<= numSites)
        allowedSz = totSz
    else
        @assert all(-numSites <= totSz <= numSites)
        allowedSz = [totSz]
    end

    # create dictionary to store configurations classified by their total occupancies.
    # For eg, basisStates = {0: [[0, 0]], 1: [[1, 0], [0, 1]], 2: [[1, 1]]}}
    basisStates = Dict{Tuple{Int64,Int64},Vector{Dict{BitVector, Float64}}}()

    # generate all possible configs, by running integers from 0 to 2^N-1, and converting
    # them to binary forms
    allConfigs = digits.(0:2^numLevels-1, base=2, pad=numLevels) |> reverse

    # loop over allowed occupancies, and store the configs in the dict
    # classified by their total occupancies. We will also store the subdimensions
    # while we are creating the dictionary.

    for (totSz, totOcc) in Iterators.product(allowedSz, allowedOccupancies)
        totOccArr = sum.(allConfigs)
        totSzArr = [sum(config[1:2:end]) - sum(config[2:2:end]) for config in allConfigs]
        allowedConfigs = allConfigs[(totOccArr .== totOcc) .& (totSzArr .== totSz)]
        for config in allowedConfigs
            if localOccupancy[2] != -1 && sum(config[localOccupancy[1]]) ≠ localOccupancy[2]
                continue
            end
            sigmaz = sum(config[1:2:end]) - sum(config[2:2:end])
            if !((totOcc, sigmaz) in keys(basisStates))
                basisStates[(totOcc, sigmaz)] = []
            end
            if length(config) > 0
                push!(basisStates[(totOcc, sigmaz)], Dict(config => 1.0))
            end
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


function operatorMatrix(basisStates::Dict{Tuple{Int64,Int64}, Vector{Dict{BitVector, Float64}}}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    operatorFullMatrix = Dict(key => zeros(length(value), length(value)) for (key, value) in basisStates)
    Threads.@threads for (key, bstates) in collect(basisStates)
        Threads.@threads for (col, colState) in collect(enumerate(bstates))
            newStateDict = applyOperatorOnState(colState, operatorList)
            for (row, rowState) in enumerate(bstates)
                if operatorFullMatrix[key][row, col] != 0
                    continue
                end
                overlap = sum([newStateDict[key] * rowState[key] for key in intersect(keys(newStateDict), keys(rowState))])
                operatorFullMatrix[key][row, col] = overlap
                operatorFullMatrix[key][col, row] = overlap'
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
