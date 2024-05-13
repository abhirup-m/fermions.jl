using Combinatorics
using ProgressMeter
using LinearAlgebra

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
function BasisStates(numLevels::Int64; totOccupancy::Union{Vector{Int64},Nothing}=nothing)

    # if totOccupancy is not set, all occupancies from 0 to N are allowed,
    # otherwise use only the provided totOccupancy
    if isnothing(totOccupancy)
        allowedOccupancies = 0:numLevels
    else
        @assert all(0 .<= totOccupancy .<= numLevels)
        allowedOccupancies = totOccupancy
    end

    # create dictionary to store configurations classified by their total occupancies.
    # For eg, basisStates = {0: [[0, 0]], 1: [[1, 0], [0, 1]], 2: [[1, 1]]}}
    basisStates = Dict{Tuple{Int64,Int64},Vector{BitArray}}()

    # generate all possible configs, by running integers from 0 to 2^N-1, and converting
    # them to binary forms
    allConfigs = digits.(0:2^numLevels-1, base=2, pad=numLevels) |> reverse

    # loop over allowed occupancies, and store the configs in the dict
    # classified by their total occupancies. We will also store the subdimensions
    # while we are creating the dictionary.

    for totOcc in allowedOccupancies
        totoccMatchingConfigs = allConfigs[sum.(allConfigs).==totOcc]
        for config in totoccMatchingConfigs
            sigmaz = sum(config[1:2:end]) - sum(config[2:2:end])
            if !((totOcc, sigmaz) in keys(basisStates))
                basisStates[(totOcc, sigmaz)] = BitArray[]
            end
            push!(basisStates[(totOcc, sigmaz)], config)
        end
    end

    return basisStates
end

function TransformBit(qubit::Bool, operator::Char)
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

function applyOperatorOnState(stateDict::Dict{BitVector,Float64}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    # define a dictionary for the final state obtained after applying 
    # the operators on all basis states
    completeOutputState = Dict{BitVector,Float64}()
    # loop over all operator tuples within operatorList
    for ((opType, opMembers), opStrength) in pairs(operatorList)

        for (state, coefficient) in stateDict

            newCoefficient = coefficient
            newState = copy(state)

            # for each basis state, obtain a modified state after applying the operator tuple
            for (siteIndex, operator) in zip(reverse(opMembers), reverse(opType))
                newQubit, factor = TransformBit(newState[siteIndex], operator)
                if factor == 0
                    newCoefficient = 0
                    break
                end
                # calculate the fermionic exchange sign by counting the number of
                # occupied states the operator has to "hop" over
                exchangeSign = operator in ['+', '-'] ? (-1)^sum(newState[1:siteIndex-1]) : 1
                @inbounds newState[siteIndex] = newQubit
                newCoefficient *= exchangeSign * factor
            end
            if newCoefficient != 0
                if newState in keys(completeOutputState)
                    completeOutputState[newState] += opStrength * newCoefficient
                else
                    completeOutputState[newState] = opStrength * newCoefficient
                end
            end
        end
    end
    return completeOutputState
end


function generalOperatorMatrix(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    operatorFullMatrix = Dict(key => zeros(length(value), length(value)) for (key, value) in basisStates)
    for (key, bstates) in collect(basisStates)
        for (index, state) in collect(enumerate(bstates))
            for (nState, coeff) in applyOperatorOnState(Dict(state => 1.0), operatorList)
                operatorFullMatrix[key][index, bstates.==[nState]] .= coeff
            end
        end
    end
    return operatorFullMatrix
end


function broadcastCouplingSet(matrixSet::Vector{Dict{Tuple{Int64,Int64}, Matrix{Float64}}}, couplingSet::Vector{Float64})
    matrixSetCouplingMultiplied = fetch.([Threads.@spawn ((matrix, coupling) -> Dict(k => v .* coupling for (k, v) in matrix))(matrix, coupling) for (matrix, coupling) in zip(matrixSet, couplingSet)])
    return merge(+, matrixSetCouplingMultiplied...)
end


function generalOperatorMatrix(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, operatorList::Vector{Tuple{String,Vector{Int64}}}, couplingMatrix::Vector{Vector{Float64}})

    # calculates the matrix form for the skeleton operator `key`. Itself is
    # a dict, with keys representing quantum numbers (ntot, Sztot) and values
    # representing matrices for the subspace corresponding to the quantum numbers.
    matrixSet = fetch.([Threads.@spawn generalOperatorMatrix(basisStates, Dict(operator => 1.0)) for operator in operatorList])
    operatorMatrixSet = fetch.([Threads.@spawn broadcastCouplingSet(matrixSet, couplingSet) for couplingSet in couplingMatrix])
    return operatorMatrixSet

end


function getSpectrum(hamiltonian::Dict{Tuple{Int64,Int64},Matrix{Float64}}; tolerance=1e-10, progressEnabled=false)
    eigvals = Dict{Tuple{Int64,Int64},Vector{Float64}}(index => [] for index in keys(hamiltonian))
    eigvecs = Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}(index => [] for index in keys(hamiltonian))
    pbar = nothing
    if progressEnabled
        pbar = Progress(length(hamiltonian))
    end
    for (index, matrix) in collect(hamiltonian)
        roundedMatrix = round.(matrix, digits=trunc(Int, -log10(tolerance)))
        F = eigen(Hermitian(roundedMatrix))
        eigvalues = F.values
        eigvals[index] = sort(eigvalues)
        eigvecs[index] = [F.vectors[:, i] for i in sortperm(eigvalues)]
        if !isnothing(pbar)
            next!(pbar)
        end
    end
    if !isnothing(pbar)
        finish!(pbar)
    end
    return eigvals, eigvecs
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
