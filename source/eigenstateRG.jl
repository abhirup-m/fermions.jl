include("./fermionise.jl")


function InitWavefunction(basisStates, hamiltonian, displayGstate=false)
    eigvals, eigstates = diagonalise(basisStates, hamiltonian)
    println(eigvals[1:2])
    # if displayGstate:
    #     print (visualise_state(basisStates, eigstates[1]))
    return eigstates[1]
end


function ApplyInverseTransform(oldState, numEntangled, unitaryOperatorFunction, alpha, sectors)
    @assert sectors in ["p", "h", "ph", "hp"] "Provided IOM sectors not among p, h or ph/hp."
    if sectors == "p"
        oldStateEnlarged = Dict([(key * "11", val) for (key, val) in oldState])
    elseif sectors == "h"
        oldStateEnlarged = Dict([(key * "00", val) for (key, val) in oldState])
    elseif sectors == "ph"
        oldStateEnlarged = Dict([(key * "1100", val) for (key, val) in oldState])
    else
        oldStateEnlarged = Dict([(key * "0011", val) for (key, val) in oldState])
    end
    unitaryOperatorList = unitaryOperatorFunction(alpha, numEntangled, sectors)
    newState = copy(oldStateEnlarged)
    batchSize = 100
    @showprogress for (basisState, coeff) in oldStateEnlarged
        for operatorBatch in Iterators.partition(unitaryOperatorList, batchSize)
            newState = merge(+, newState, merge(+, pmap(ApplyOperatorOnState, [(Dict([(basisState, coeff)]), [operator])
                                                                               for operator in operatorBatch])...))
        end
    end

    total_norm = sum(values(newState) .^ 2)^0.5
    map!(x -> x / total_norm, values(newState))

    return newState
end

function getWavefunctionRG(initCouplings, alphaArray, numEntangled, numReverseSteps, hamiltonianFunction, unitaryOperatorFunction, sectors, displayGstate=false)

    @assert length(alphaArray) >= numReverseSteps "Number of values of 'alpha' is not enough for the requested number of reverse RG steps."

    basisStates = BasisStates(2 * (1 + numEntangled))
    hamiltonianMatrix = hamiltonianFunction(basisStates, numEntangled, initCouplings)
    initState = InitWavefunction(basisStates, hamiltonianMatrix)
    stateFlowArray = [initState]

    for (i, alpha) in enumerate(alphaArray[1:numReverseSteps])
        newState = ApplyInverseTransform(stateFlowArray[end], numEntangled + length(sectors) * (i - 1), unitaryOperatorFunction, alpha, sectors)
        push!(stateFlowArray, newState)
    end

    return stateFlowArray
end
