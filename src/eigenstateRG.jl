function getWavefunctionRG(initState::Dict{BitVector, Float64}, alphaValues::Vector{Float64}, numSteps::Integer, unitaryOperatorFunction::Function, stateExpansionFunction::Function, sectors::String; tolerance::Float64=1e-16)
    
    numEntangled = div(length(collect(keys(initState))[1]), 2)
    stateFlowArray = Dict{BitVector, Float64}[]
    push!(stateFlowArray, initState)
    
    pbar = Progress(numSteps; dt=0.1)
    for alpha in alphaValues
        newState = stateExpansionFunction(stateFlowArray[end])
        unitaryOperatorList = unitaryOperatorFunction(alpha, numEntangled, sectors)
        numEntangled = div(length(collect(keys(newState))[1]), 2)
        stateRenormalisation = fetch.([Threads.@spawn ApplyOperator(newState, Dict(k => v)) for (k,v) in unitaryOperatorList])
        mergewith!(+, newState, stateRenormalisation...)

        newState= Dict(k => v for (k,v) in newState if abs(v) > tolerance)
        total_norm = sum(values(newState) .^ 2)^0.5
        newState= Dict(k => v/total_norm for (k,v) in newState)
        push!(stateFlowArray, newState)
        next!(pbar)
    end

    return stateFlowArray
end
