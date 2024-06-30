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
function BasisStates(numLevels::Int64; allowedSz::Union{Int64,Vector{Int64},Nothing}=nothing, allowedOccupancies::Union{Float64,Vector{Float64},Nothing}=nothing, localOccupancy::Tuple{Vector{Int64}, Int64}=([1,2], -1))
    @assert numLevels > 0 && numLevels % 2 == 0

    # if totOccupancy is not set, all occupancies from 0 to N are allowed,
    # otherwise use only the provided totOccupancy
    if !isnothing(allowedOccupancies)
        if typeof(allowedOccupancies) == Float64
            allowedOccupancies = [allowedOccupancies]
        end
        @assert all(0 .<= allowedOccupancies .<= numLevels)
    end

    # same for total Sz
    numSites = trunc(Int, numLevels / 2)
    if !isnothing(allowedSz)
        if typeof(allowedSz) == Int64
            allowedSz = [allowedSz]
        end
        @assert all(-numSites .<= allowedSz .<= numSites)
    end

    # create dictionary to store configurations classified by their total occupancies.
    # For eg, basisStates = {0: [[0, 0]], 1: [[1, 0], [0, 1]], 2: [[1, 1]]}}
    basisStates = Dict{Tuple{Float64,Int64},Vector{Dict{BitVector, Float64}}}()

    # generate all possible configs, by running integers from 0 to 2^N-1, and converting
    # them to binary forms
    allConfigs = digits.(0:2^numLevels-1, base=2, pad=numLevels) |> reverse

    for config in allConfigs

        # check which states have the correct occupancy and total Sz
        totOcc = sum(config) / numLevels
        totSz = sum(config[1:2:end]) - sum(config[2:2:end])
        if !isnothing(allowedOccupancies) && totOcc ∉ allowedOccupancies
            continue
        end
        if !isnothing(allowedSz) && totSz ∉ allowedSz
            continue
        end

        # for the ones that do, check if any local occupancy criteria is satisfied
        if localOccupancy[2] != -1 && sum(config[localOccupancy[1]]) ≠ localOccupancy[2]
            continue
        end

        if !((totOcc, totSz) in keys(basisStates))
            basisStates[(totOcc, totSz)] = []
        end
        if length(config) > 0
            push!(basisStates[(totOcc, totSz)], Dict(config => 1.0))
        end
    end

    return basisStates
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


function applyOperatorOnState(stateDict::Dict{BitVector,Float64}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    @assert maximum([maximum(opMembers) for (_, opMembers) in keys(operatorList)]) ≤ length(collect(keys(stateDict))[1])

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
                    @inbounds completeOutputState[newState] += opStrength * newCoefficient
                else
                    @inbounds completeOutputState[newState] = opStrength * newCoefficient
                end
            end
        end
    end
    return completeOutputState
end


"""Creates matrix representation for the provided abstract operator definition.
"""
function operatorMatrix(basisStates::Dict{Tuple{Float64,Int64}, Vector{Dict{BitVector, Float64}}}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    operatorFullMatrix = Dict(key => zeros(length(value), length(value)) for (key, value) in basisStates)

    # loop over symmetry sectors of the basis
    Threads.@threads for (key, bstates) in collect(basisStates)

        # loop over basis states within the symmetry sector.
        # These states will act as incoming states in terms
        # of the operator (hence the column index 'col').
        Threads.@threads for (col, colState) in collect(enumerate(bstates))

            # get action of operator on incoming basis state |ψ>
            newStateDict = applyOperatorOnState(colState, operatorList)

            # again loop over the basis states, this term for the outgoing states <ϕ|
            for (row, rowState) in enumerate(bstates)
                if operatorFullMatrix[key][row, col] != 0
                    continue
                end

                # calculate overlap between the outgoing state <ϕ| and the 
                # modified incoming state O|ψ>.
                overlap = sum([newStateDict[key] * rowState[key] for key in intersect(keys(newStateDict), keys(rowState))])
                operatorFullMatrix[key][row, col] = overlap
                operatorFullMatrix[key][col, row] = overlap'
            end
        end
    end
    return operatorFullMatrix
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


function getSpectrum(
        basisStates::Dict{Tuple{Float64,Int64},Vector{Dict{BitVector, Float64}}},
        hamiltonian::Dict{Tuple{Float64,Int64},Matrix{Float64}}; 
        tolerance::Float64=1e-16, maxNum::Int64=0, progressEnabled::Bool=false
    )
    @assert all(x -> ishermitian(x), values(hamiltonian))
    eigvals = Dict{Tuple{Float64,Int64},Vector{Float64}}(index => [] for index in keys(hamiltonian))
    eigvecs = Dict{Tuple{Float64,Int64},Vector{Dict{BitVector, Float64}}}(index => [] for index in keys(hamiltonian))
    pbar = nothing
    if progressEnabled
        pbar = Progress(length(hamiltonian))
    end
    spectra = Dict(zip(keys(hamiltonian), 
                       fetch.([Threads.@spawn eigen(Hermitian(round.(matrix, digits=trunc(Int, -log10(tolerance)))))
                               for matrix in values(hamiltonian)])
                      )
                  )
    for (index, F) in spectra
        eigvalues = F.values
        maxNum = ifelse(maxNum == 0, length(eigvalues), minimum((length(eigvalues), maxNum)))
        @inbounds eigvals[index] = sort(eigvalues)[1:maxNum]
        @inbounds eigvecs[index] = rotateStates(basisStates[index], 
                                                [round.(F.vectors[:, i], digits=trunc(Int, -log10(tolerance))) 
                                                 for i in sortperm(eigvalues)[1:maxNum]])
        
        if !isnothing(pbar)
            next!(pbar)
        end
    end
    if !isnothing(pbar)
        finish!(pbar)
    end
    return eigvals, eigvecs
end


function rotateStates(basis::Vector{Dict{BitVector, Float64}}, transformation::Vector{Vector{Float64}})
    newBasis = Dict{BitVector, Float64}[Dict{BitVector, Float64}() for _ in transformation]
    # loop over the eigenvectors of the transformation
    # (for eg. eigenvector = [-0.5, 0.5, 0.5, 0.5])
    for (i, eigenvector) in collect(enumerate(transformation))

        # loop over the v_i^j
        for (j, multiplier) in enumerate(eigenvector)
            if abs(multiplier) < 1e-10
                continue
            end

            # create the new basis states: {|1>: v_i^1 * b_{1, old}, |2>: v_i^2 * b_{2, old} ... }
            multipliedState = Dict(q => multiplier * b for (q, b) in basis[j])
            mergewith!(+, newBasis[i], multipliedState)
        end
        newBasis[i] = Dict(k => v for (k,v) in newBasis[i] if abs(v) > 1e-10)
        @assert !isempty(newBasis[i])
    end
    return newBasis
end
