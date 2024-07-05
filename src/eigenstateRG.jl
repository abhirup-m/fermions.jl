function getWavefunctionRG(initState::Dict{BitVector, Float64}, numEntangledSites::Integer, numReverseSteps::Integer, unitaryOperatorFunction, stateExpansionFunction, sectors::String; tolerance::Float64=1e-16)
    
    stateFlowArray = Dict{BitVector, Float64}[]
    push!(stateFlowArray, initState)
    
    pbar = Progress(numReverseSteps; dt=0.1)
    for step in 1:numReverseSteps
        newState = stateExpansionFunction(stateFlowArray[end])
        unitaryOperatorList = unitaryOperatorFunction(step)
        stateRenormalisation = fetch.([Threads.@spawn applyOperatorOnState(newState, Dict(k => v)) for (k,v) in unitaryOperatorList])
        mergewith!(+, newState, stateRenormalisation...)

        newState= Dict(k => v for (k,v) in newState if abs(v) > tolerance)
        total_norm = sum(values(newState) .^ 2)^0.5
        newState= Dict(k => v/total_norm for (k,v) in newState)
        push!(stateFlowArray, newState)
        next!(pbar)
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
