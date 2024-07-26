using Combinatorics

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


# rounds off a given float to the desired tolerance
function roundTo(val, tolerance)
    return round(val, digits=trunc(Int, -log10(tolerance)))
end


function BasisStates(numLevels::Int64, totOcc::Int64; magzCriteria::Function=x -> true, localCriteria::Function=x -> true)
    basis = Dict{BitVector,Float64}[]
    configs = Combinatorics.multiset_permutations([fill(1, totOcc); fill(0, numLevels - totOcc)], numLevels) |> collect
    for config in configs
        if magzCriteria(config) && localCriteria(config)
            push!(basis, Dict(BitVector(config) => 1.0))
        end
    end
    return basis
end


function BasisStates(numLevels::Int64; totOccCriteria::Function=(x, N) -> true, magzCriteria::Function=x -> true, localCriteria::Function=x -> true)
    basis = Dict{BitVector,Float64}[]
    for decimalNum in 0:2^numLevels-1
        config = digits(decimalNum, base=2, pad=numLevels) |> reverse
        if totOccCriteria(config, numLevels) && magzCriteria(config) && localCriteria(config)
            push!(basis, Dict(BitVector(config) => 1.0))
        end
    end
    return basis
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


function ApplyOperator(operator::Vector{Tuple{String,Vector{Int64},Float64}}, incomingState::Dict{BitVector,Float64}; addOn=false, tolerance=1e-16)
    @assert !isempty(operator)
    @assert maximum([maximum(positions) for (_, positions, _) in operator]) ≤ length.(keys(incomingState))[1]

    replaceState = ifelse(addOn, copy(incomingState), Dict{BitVector,Float64}())
    if addOn
        map!(x -> x/length(operator), values(replaceState))
    end
    outgoingState = [copy(replaceState) for i in 1:length(operator)]
    # loop over all operator tuples within operatorList
    Threads.@threads for i in eachindex(operator)
        goodIncomingStates = copy(incomingState)
        opType, opMembers, opStrength = operator[i]
        for i in eachindex(opMembers)[end:-1:2]
            if opMembers[i] ∈ opMembers[i+1:end]
                continue
            end
            if opType[i] == '+' || opType[i] == 'h'
                filter!(p -> p[1][opMembers[i]] == 0, goodIncomingStates)
            else
                filter!(p -> p[1][opMembers[i]] == 1, goodIncomingStates)
            end
        end
        for (incomingBasisState, coefficient) in goodIncomingStates

            newCoefficient = coefficient
            outgoingBasisState = copy(incomingBasisState)

            # for each basis state, obtain a modified state after applying the operator tuple
            for (siteIndex, operator) in zip(reverse(opMembers), reverse(opType))
                newQubit, factor = TransformBit(outgoingBasisState[siteIndex], operator)
                if factor == 0
                    newCoefficient = 0
                    break
                end
                # calculate the fermionic exchange sign by counting the number of
                # occupied states the operator has to "hop" over
                exchangeSign = ifelse(operator in ['+', '-'], (-1)^sum(outgoingBasisState[1:siteIndex-1]), 1)

                outgoingBasisState[siteIndex] = newQubit
                newCoefficient *= exchangeSign * factor
            end

            if abs(newCoefficient) > tolerance
                if haskey(outgoingState[i], outgoingBasisState)
                    outgoingState[i][outgoingBasisState] += opStrength * newCoefficient
                else
                    outgoingState[i][outgoingBasisState] = opStrength * newCoefficient
                end
            end
        end
    end

    return mergewith(+, outgoingState...)
end


function OperatorMatrix(basisStates::Vector{Dict{BitVector,Float64}}, operator::Vector{Tuple{String,Vector{Int64},Float64}})
    operatorMatrix = zeros(length(basisStates), length(basisStates))
    newStates = fetch.([Threads.@spawn ApplyOperator(operator, incomingState) for incomingState in basisStates])
    Threads.@threads for incomingIndex in eachindex(newStates)
        Threads.@threads for outgoingIndex in eachindex(basisStates)
            operatorMatrix[outgoingIndex, incomingIndex] = StateOverlap(basisStates[outgoingIndex], newStates[incomingIndex])
        end
    end
    return operatorMatrix
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


function StateOverlap(state1::Dict{BitVector,Float64}, state2::Dict{BitVector,Float64})
    overlap = 0
    keys2 = keys(state2)
    for (key, val) in state1
        if key ∈ keys2
            overlap += val * state2[key]
        end
    end
    return overlap
end
