function getWavefunctionRG(initState::Dict{BitVector, Float64}, initCouplings, numEntangledSites::Integer, numReverseSteps::Integer, hamiltonianFunction, unitaryOperatorFunction, stateExpansionFunction, sectors::String, displayGstate=false)
    
    stateFlowArray = Dict{BitVector, Float64}[]
    push!(stateFlowArray, initState)

    for step in 1:numReverseSteps
        newState = stateExpansionFunction(stateFlowArray[end])
        unitaryOperatorList = unitaryOperatorFunction(step)
        stateRenormalisation = fetch.([Threads.@spawn applyOperatorOnState(newState, Dict(k => v)) for (k,v) in unitaryOperatorList])
        mergewith!(+, newState, stateRenormalisation...)

        total_norm = sum(values(newState) .^ 2)^0.5
        map!(x -> x / total_norm, values(newState))
        push!(stateFlowArray, newState)
    end

    return stateFlowArray
end


function correlationRG(stateFlowArray::Vector{Dict{BitVector, Float64}}, correlationOperators)
    numLevelsArr = [length(collect(keys(state))[1]) for state in stateFlowArray]
    correlationRGResults = [Float64[] for _ in correlationOperators]
    for state in stateFlowArray
        for (i, operator) in enumerate(correlationOperators)
            push!(correlationRGResults[i], 0)
            MPsi = applyOperatorOnState(state, operator)
            for (bstate, coeff) in MPsi
                correlationRGResults[i][end] += bstate ∈ keys(state) ? coeff * state[bstate] : 0
            end
        end
    end
    return correlationRGResults
end
