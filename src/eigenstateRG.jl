function InitWavefunction(basisStates::Vector{BitVector}, initCouplings, hamiltonianFunction; tolerance=1e-10)
    operatorList = hamiltonianFunction(initCouplings...)
    hamiltonianMatrix = generalOperatorMatrix(basisStates, operatorList)
    groundStates = getGstate(hamiltonianMatrix)
    @assert length(groundStates) == 1
    stateDict = Dict(bstate => coeff for (bstate,coeff) in zip(basisStates, groundStates[1]) if abs(coeff) > tolerance)
    return stateDict
end


function ApplyInverseTransform(state, unitaryOperatorFunction, alpha, sectors)
    @assert sectors in ["p", "h", "ph"] "Provided IOM sectors not among p, h or ph/hp."
    numEntangled = trunc(Int, length(collect(keys(state))[1]) / 2)
    if sectors == "p"
        state = Dict{keytype(state), valtype(state)}([key; [1, 1]] => val for (key, val) in state)
    elseif sectors == "h"
        state = Dict{keytype(state), valtype(state)}([key; [0, 0]] => val for (key, val) in state)
    else
        state = Dict{keytype(state), valtype(state)}([key; [1, 1, 0, 0]] => val for (key, val) in state)
    end
    unitaryOperatorList = unitaryOperatorFunction(alpha, numEntangled, sectors)
    stateRenormalisation = fetch.([Threads.@spawn applyOperatorOnState(state, Dict(k => v)) for (k,v) in unitaryOperatorList])
    mergewith!(+, state, stateRenormalisation...)

    total_norm = sum(values(state) .^ 2)^0.5
    map!(x -> x / total_norm, values(state))

    return state
end

function getWavefunctionRG(basisFunction, initCouplings, alphaArray, numEntangled::Integer, numReverseSteps::Integer, hamiltonianFunction, unitaryOperatorFunction, sectors::String, displayGstate=false)
    
    @assert length(alphaArray) ≥ numReverseSteps "Number of values of 'alpha' is not enough for the requested number of reverse RG steps."
    basisStates = basisFunction(2 * numEntangled + 2)
    @assert length(basisStates) == 1
    basisStates = collect(values(basisStates))[1]
    initState = InitWavefunction(basisStates, initCouplings, hamiltonianFunction)
    stateFlowArray = Dict{BitVector, Float64}[]
    push!(stateFlowArray, initState)

    for (i, alpha) in enumerate(alphaArray[1:numReverseSteps])
        newState = ApplyInverseTransform(stateFlowArray[end], unitaryOperatorFunction, alpha, sectors)
        push!(stateFlowArray, newState)
    end

    return stateFlowArray
end


function correlationRG(basisFunction, stateFlowArray::Vector{Dict{BitVector, Float64}}, correlationOperators)
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
