using Combinatorics
using Distributed
using ProgressMeter
using SparseArrays


function BasisStates(numLevels, totOccupancy=-999, totSpin=-999)
	basisStates = []
	for i in 0:2^(numLevels)-1
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


function TransformQubit(qubit, operator)
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


function GeneralOperatorMatrix(basisStates, operatorList, show_progress=true)
    dimension = length(basisStates)
    operatorMatrix = zeros(dimension, dimension)
    pbar = Progress(length(basisStates))
    Threads.@threads for (index, state) in collect(enumerate(basisStates))
        newState = zeros(length(state))
        for (couplingStrength, operatorChunk, siteIndices) in operatorList
            newState[:] = state
            matrixElement = 1
            for (siteIndex, operator) in zip(reverse(siteIndices), reverse(operatorChunk))
                newQubit, factor = TransformQubit(state[siteIndex], operator)
                newState[siteIndex] = newQubit
                exchangeSign = operator in ['+', '-'] ? (-1)^sum(state[1:siteIndex]) : 1
                matrixElement *= factor * exchangeSign
            end
            if matrixElement != 0
                newIndex = findfirst(x->x==newState, basisStates)
                operatorMatrix[index, newIndex] += matrixElement * couplingStrength
            end
        end
        next!(pbar)
    end
    return operatorMatrix
end    