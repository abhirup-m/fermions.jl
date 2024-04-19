using LinearAlgebra
using Distributed
addprocs(1)
@everywhere using Combinatorics
@everywhere using ProgressMeter
# @everywhere using Term.Progress
@everywhere using OrderedCollections


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
function BasisStates(numLevels::Int64; totOccupancy::Union{Vector{Int64},Nothing}=nothing)

    # if totOccupancy is not set, all occupancies from 0 to N are allowed,
    # otherwise use only the provided totOccupancy
    if isnothing(totOccupancy)
        allowedOccupancies = 0:numLevels
    else
        @assert all(0 .<= totOccupancy .<= numLevels)
        allowedOccupancies = totOccupancy
    end

    # create dictionary to store configurations classified by their total occupancies.
    # For eg, basisStates = {0: [[0, 0]], 1: [[1, 0], [0, 1]], 2: [[1, 1]]}}
    basisStates = OrderedDict((totOcc => BitArray[] for totOcc in allowedOccupancies))

    # generate all possible configs, by running integers from 0 to 2^N-1, and converting
    # them to binary forms
    allConfigs = digits.(0:2^numLevels-1, base=2, pad=numLevels) |> reverse

    # loop over allowed occupancies, and store the configs in the dict
    # classified by their total occupancies. We will also store the subdimensions
    # while we are creating the dictionary.
    subDimensions = Vector{Int64}(undef, length(allowedOccupancies))
    Threads.@threads for index in eachindex(allowedOccupancies)
        totOcc = allowedOccupancies[index]
        basisStates[totOcc] = allConfigs[sum.(allConfigs).==totOcc]
        subDimensions[index] = length(basisStates[totOcc])
    end

    return basisStates, subDimensions
end

@everywhere function TransformBit(qubit::Bool, operator::Char)
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

function applyOperatorOnState(stateDict::Dict{BitVector,Float64}, operatorList::Vector{Tuple{String,Float64,Vector{Int64}}})
    # define a dictionary for the final state obtained after applying 
    # the operators on all basis states
    completeOutputState = Dict{BitVector,Float64}()

    # loop over all operator tuples within operatorList
    for (opType, opStrength, opMembers) in operatorList
        @assert unique(opMembers) == opMembers

        holeProjectionSites = opMembers[findall(x -> x == '+' || x == 'h', opType)]
        particleProjectionSites = opMembers[findall(x -> x == '-' || x == 'n', opType)]

        for (state, coefficient) in stateDict

            # from the outset, ignore states that are creating at occupied bits, destroying 
            # at empty bits, applying n at empty bits or checking 1-n at occupied bits.
            if 1 in state[holeProjectionSites] || 0 in state[particleProjectionSites]
                continue
            end

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

            # if the coefficient of the transformation is non-zero, add this to the output state
            # dictionary for the operator tuple we are currently looping over
            if newCoefficient != 0
                if newState in keys(completeOutputState)
                    completeOutputState[newState] += opStrength * newCoefficient
                else
                    completeOutputState[newState] = opStrength * newCoefficient
                end
            end
        end
    end
    return completeOutputState
end

function generalOperatorMatrix(basisStates::OrderedDict{Int64,Vector{BitArray}}, operatorList::Vector{Tuple{String,Float64,Vector{Int64}}})
    operatorFullMatrix = OrderedDict{Int64,Matrix{Float64}}()
    Threads.@threads for (totOcc, bstates) in collect(pairs(basisStates))
        operatorSubMatrix = zeros(length(bstates), length(bstates))
        Threads.@threads for (index, state) in collect(enumerate(bstates))
            newState = applyOperatorOnState(Dict(state => 1.0), operatorList)
            matchingIndices = findall(in(keys(newState)),bstates)
            @inbounds operatorSubMatrix[index, matchingIndices] .= collect(values(newState)) 
        end
        operatorFullMatrix[totOcc] = operatorSubMatrix
    end
    return operatorFullMatrix
end
