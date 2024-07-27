"""
    roundTo(val, tolerance)

Round the provided float to the specified tolerance.
# Examples
```jldoctest
julia> fermions.roundTo(1.122323, 1e-3)
1.122
```
"""
function roundTo(
        val::Float64, 
        tolerance::Float64
    )
    return round(val, digits=trunc(Int, -log10(tolerance)))
end


function BasisStates(
        numLevels::Int64, 
        totOccReq::Vector{Int64},
        magzReq::Vector{Int64},
        localCriteria::Function
    )
    @assert !isempty(totOccReq) && !isempty(magzReq)
    basis = Dict{BitVector,Float64}[]
    for decimalNum in 0:2^numLevels-1
        config = digits(decimalNum, base=2, pad=numLevels) |> reverse
        if !isempty(totOccReq)
            totOcc = sum(config)
            if totOcc ∉ totOccReq
                continue
            end
        end
        if !isempty(magzReq)
            magz = sum(config[1:2:end]) - sum(config[2:2:end])
            if magz ∉ magzReq
                continue
            end
        end
        if localCriteria(config)
            push!(basis, Dict(BitVector(config) => 1.0))
        end
    end
    return basis
end


function BasisStates(
        numLevels::Int64;
        totOccReq::Union{Vector{Int64}, Int64}=nothing,
        magzReq::Union{Vector{Int64}, Int64}=nothing,
        localCriteria::Function=x -> true
    )
    if isnothing(totOccReq)
        totOccReq = collect(0, numLevels)
    elseif typeof(totOccReq) == Int64
        totOccReq = [totOccReq]
    end
    if isnothing(magzReq)
        magzReq = collect(numLevels - div(numLevels, 2), -div(numLevels, 2))
    elseif typeof(magzReq) == Int64
        magzReq = [magzReq]
    end
    return BasisStates(numLevels, totOccReq, magzReq, localCriteria)
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
