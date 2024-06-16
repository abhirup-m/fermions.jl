function InitWavefunction(basisStates, initCouplings, hamiltonianFunction)
    operatorList = hamiltonianFunction(initCouplings...)
    hamiltonianMatrix = generalOperatorMatrix(basisStates, operatorList)
    eigvals, eigstates = getSpectrum(hamiltonianMatrix)
    gsEnergy, groundStates, blocks = getGstate(eigvals, eigstates)
    @assert length(blocks) == 1
    return Dict(zip(basisStates[blocks[1]], groundStates[1]))
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

function getWavefunctionRG(initCouplings, alphaArray, numEntangled::Integer, numReverseSteps::Integer, hamiltonianFunction, unitaryOperatorFunction, sectors::String, displayGstate=false)
    
    @assert length(alphaArray) ≥ numReverseSteps "Number of values of 'alpha' is not enough for the requested number of reverse RG steps."

    basisStates = BasisStates(2 * (1 + numEntangled))
    initState = InitWavefunction(basisStates, initCouplings, hamiltonianFunction)
    stateFlowArray = [initState]

    for (i, alpha) in enumerate(alphaArray[1:numReverseSteps])
        newState = ApplyInverseTransform(stateFlowArray[end], unitaryOperatorFunction, alpha, sectors)
        push!(stateFlowArray, newState)
    end

    return stateFlowArray
end
