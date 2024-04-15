using Combinatorics
using Distributed
using ProgressMeter
using SparseArrays


function BasisStates(numLevels::Int64; totOccupancy::Float64=-999., totSpin::Float64=-999.)
    basisStates = Vector{Int64}[]
	Threads.@threads for i in 0:2^(numLevels)-1
		state = [parse(Int, ch) for ch in lpad(string(i, base=2), numLevels, '0')]
        stateTotOccupancy = sum(state)
        stateTotSpin = 0.5 * (sum(state[1:2:end]) - sum(state[2:2:end]))
		if ((totOccupancy == -999 || stateTotOccupancy == totOccupancy)
            && (totSpin == -999 || stateTotSpin == totSpin))
            push!(basisStates, state)
		end
	end
	return basisStates
end


function TransformQubit(qubit::Int64, operator::Char)
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


function applyOperatorOnState(stateDict::Dict{Vector{Int64}, Float64}, operatorList::Vector{Tuple{String, Float64, Vector{Int64}}})
    completeOutputState = @distributed (a, b) -> merge(+, a, b) for (opType, opStrength, opMembers) in operatorList
        outputState = Dict()
        for (state, coefficient) in stateDict
            newCoefficient = coefficient
            newState = copy(state)
            for (siteIndex, operator) in zip(reverse(opMembers), reverse(opType))
                newQubit, factor = TransformQubit(newState[siteIndex], operator)
                if factor == 0
                    newCoefficient = 0
                    break
                end
                newState[siteIndex] = newQubit
                exchangeSign = operator in ['+', '-'] ? (-1)^sum(newState[1:siteIndex]) : 1
                newCoefficient *= factor
            end
            if newCoefficient != 0
                if newState in keys(outputState)
                    outputState[newState] += newCoefficient
                else
                    outputState[newState] = newCoefficient
                end
            end
        end
        outputState;
    end
    return completeOutputState
end

function generalOperatorMatrix(basisStates::Vector{Vector{Int64}}, operatorList::Vector{Tuple{String, Float64, Vector{Int64}}})
    operatorDimension = length(basisStates)
    operatorMatrix = zeros((operatorDimension, operatorDimension))
    #=Threads.@threads =#
    for (index1, state) in collect(enumerate(basisStates))
        stateDict = Dict(state => 1.)
        for (newState, coefficient) in applyOperatorOnState(stateDict, operatorList)
            operatorMatrix[index1, basisStates .== [newState]] .= coefficient
        end
    end
    return operatorMatrix
end
