using Combinatorics
using Distributed
using ProgressMeter
using SparseArrays


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
function BasisStates(numLevels::Int64; totOccupancy::Int64=-999, totSpin::Float64=-999.0)

    # if totOccupancy is not set, all occupancies from 0 to N are allowed,
    # otherwise use only the provided totOccupancy
    if totOccupancy == -999
        allowedOccupancies = 0:numLevels
    else
        @assert 0 <= totOccupancy <= numLevels
        allowedOccupancies = (totOccupancy,)
    end

    # create dictionary to store configurations classified by their total occupancies.
    # For eg, basisStates = {0: [[0, 0]], 1: [[1, 0], [0, 1]], 2: [[1, 1]]}}
    basisStates = Dict((totOcc => BitArray[] for totOcc in allowedOccupancies))

    # generate all possible configs, by running integers from 0 to 2^N-1, and converting
    # them to binary forms
    allConfigs = digits.(0:2^numLevels - 1, base=2, pad=numLevels) |> reverse

    # loop over allowed occupancies, and store the configs in the dict
    # classified by their total occupancies.
    Threads.@threads for totOcc in allowedOccupancies
        basisStates[totOcc] = allConfigs[sum.(allConfigs) .== totOcc]
    end

    return basisStates
end


"""
Transforms a single Fock state |n> (n=0/1) given an operator (= c/c^†/n̂/1-n̂).
"""
function TransformBit(qubit::Int64, operator::Char)
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


function applyOperatorOnState(stateDict::Dict{Vector{Int64},Float64}, operatorList::Vector{Tuple{String,Float64,Vector{Int64}}})
    completeOutputState = @distributed (a, b) -> merge(+, a, b) for (opType, opStrength, opMembers) in operatorList
        outputState = Dict()
        for (state, coefficient) in stateDict
            newCoefficient = coefficient
            newState = copy(state)
            for (siteIndex, operator) in zip(reverse(opMembers), reverse(opType))
                newQubit, factor = TransformBit(newState[siteIndex], operator)
                if factor == 0
                    newCoefficient = 0
                    break
                end
                newState[siteIndex] = newQubit
                exchangeSign = operator in ['+', '-'] ? (-1)^sum(newState[1:siteIndex]) : 1
                newCoefficient *= exchangeSign * factor
            end
            if newCoefficient != 0
                if newState in keys(outputState)
                    outputState[newState] += opStrength * newCoefficient
                else
                    outputState[newState] = opStrength * newCoefficient
                end
            end
        end
        outputState
    end
    return completeOutputState
end

function generalOperatorMatrix(basisStates::Vector{Vector{Int64}}, operatorList::Vector{Tuple{String,Float64,Vector{Int64}}})
    operatorDimension = length(basisStates)
    operatorMatrix = zeros((operatorDimension, operatorDimension))
    Threads.@threads for (index1, state) in collect(enumerate(basisStates))
        stateDict = Dict(state => 1.0)
        for (newState, coefficient) in applyOperatorOnState(stateDict, operatorList)
            operatorMatrix[index1, basisStates.==[newState]] .= coefficient
        end
    end
    return operatorMatrix
end
