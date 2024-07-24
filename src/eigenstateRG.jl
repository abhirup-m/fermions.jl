function getWavefunctionRG(initState::Dict{BitVector,Float64}, alphaValues::Vector{Float64}, numSteps::Integer, unitaryOperatorFunction::Function, stateExpansionFunction::Function, sectors::String; tolerance::Float64=1e-16)

    @assert numSteps ≤ length(alphaValues)
    numEntangled = div(length(collect(keys(initState))[1]), 2)
    stateFlowArray = Dict{BitVector,Float64}[]
    push!(stateFlowArray, initState)

    @showprogress for alpha in alphaValues[1:numSteps]
        newState = stateExpansionFunction(stateFlowArray[end], sectors)
        unitaryOperatorList = unitaryOperatorFunction(alpha, numEntangled, sectors)
        numEntangled = div(length(collect(keys(newState))[1]), 2)
        mergewith!(+, newState, fetch.([Threads.@spawn ApplyOperator([operator], newState) for operator in unitaryOperatorList])...)


        keepIndices = abs.(values(newState)) ./ maximum(abs.(values(newState))) .> tolerance
        newState = Dict(collect(keys(newState))[keepIndices] .=> collect(values(newState))[keepIndices])
        total_norm = sum(values(newState) .^ 2)^0.5
        newState = Dict(keys(newState) .=> values(newState) ./ total_norm)
        push!(stateFlowArray, newState)
    end

    return stateFlowArray
end
