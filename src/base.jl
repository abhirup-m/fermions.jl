"""
    BasisStates(numLevels, totOccReq, magzReq, localCriteria)

Creates a set of basis states for the `numLevels` Fock states, having the
required total occupancy and magnetization, and subject to any criteria
imposed locally.

# Examples
```jldoctest
julia> BasisStates(2, [1], [0], x->true)
Dict{BitVector, Float64}[]
julia> BasisStates(3, [1], [-1], x->true)
1-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 1, 0] => 1.0)
julia> BasisStates(4, [2], [0], x->true)
4-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 0, 1, 1] => 1.0)
 Dict([0, 1, 1, 0] => 1.0)
 Dict([1, 0, 0, 1] => 1.0)
 Dict([1, 1, 0, 0] => 1.0)
julia> BasisStates(4, [2], [0], x->sum(x[1:2])==1)
2-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 1, 1, 0] => 1.0)
 Dict([1, 0, 0, 1] => 1.0)
```
"""
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
export BasisStates


"""
    BasisStates(numLevels; totOccReq, magzReq, localCriteria)

Extends BasisStates() by making the last three arguments optional.
When skipped, these arguments are filled with all possible values
and then passed to BasisStates().

# Examples
```jldoctest
julia> BasisStates(2)
4-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 0] => 1.0)
 Dict([0, 1] => 1.0)
 Dict([1, 0] => 1.0)
 Dict([1, 1] => 1.0)

julia> BasisStates(2; totOccReq=1)
2-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 1] => 1.0)
 Dict([1, 0] => 1.0)

julia> BasisStates(2; magzReq=0)
2-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 0] => 1.0)
 Dict([1, 1] => 1.0)
```
"""
function BasisStates(
        numLevels::Int64;
        totOccReq::Union{Vector{Int64}, Int64, Nothing}=nothing,
        magzReq::Union{Vector{Int64}, Int64, Nothing}=nothing,
        localCriteria::Function=x -> true
    )
    if isnothing(totOccReq)
        totOccReq = collect(0:numLevels)
    elseif typeof(totOccReq) == Int64
        totOccReq = [totOccReq]
    end
    if isnothing(magzReq)
        magzReq = collect(-div(numLevels, 2):numLevels - div(numLevels, 2))
    elseif typeof(magzReq) == Int64
        magzReq = [magzReq]
    end
    return BasisStates(numLevels, totOccReq, magzReq, localCriteria)
end
export BasisStates


function BasisStates1p(
        numLevels::Int64, 
    )
    config = BitVector(fill(0, numLevels))
    config[1] = 1
    basis = Dict{BitVector,Float64}[]
    for shift in 1:numLevels
        push!(basis, Dict(copy(config) => 1.0))
        circshift!(config, 1)
    end
    return basis
end
export BasisStates1p


"""
    TransformBit(qubit, operator)

Apply the single qubit operator ('n', 'h', '+' or '-') on a single fock state.

# Examples
```jldoctest
julia> TransformBit(Bool(0), '+')
(1, 1)

julia> TransformBit(Bool(1), 'h')
(1, 0)

julia> TransformBit(Bool(1), '-')
(0, 1)
```
"""
function TransformBit(qubit::Bool, operator::Char)
    @assert operator in ('n', 'h', '+', '-')
    if operator == 'n'
        return 0 + qubit, 0 + qubit
    elseif operator == 'h'
        return 0 + qubit, 1 - qubit
    elseif (operator == '+' && qubit == 0) || (operator == '-' && qubit == 1)
        return 1 - qubit, 1
    else
        return 0 + qubit, 0
    end
end
export TransformBit


"""
    ApplyOperatorChunk(opType, opMembers, opStrength, incomingState; tolerance)

Apply a single tensor product operator chunk (for eg., c^†_1 c_2 or n_1 c^†_3 c_4) 
on a general state and return the new state.

# Examples
```jldoctest
julia> state = Dict(Bool.([1, 0]) => 1.0, Bool.([0, 1]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 0] => 1.0
  [0, 1] => -0.5

julia> opType, opMembers, opStrength = ("+-", [1, 2], 0.1)
("+-", [1, 2], 0.1)

julia> ApplyOperatorChunk(opType, opMembers, opStrength, state)
Dict{BitVector, Float64} with 1 entry:
  [1, 0] => -0.05
```
"""
function ApplyOperatorChunk(
        opType::String,
        opMembers::Vector{Int64},
        opStrength::Float64,
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )
    outgoingState = Dict{BitVector,Float64}()
    for i in eachindex(opMembers)[end:-1:2]
        if opMembers[i] ∈ opMembers[i+1:end]
            continue
        end
        if opType[i] == '+' || opType[i] == 'h'
            filter!(p -> p[1][opMembers[i]] == 0, incomingState)
        else
            filter!(p -> p[1][opMembers[i]] == 1, incomingState)
        end
    end
    for (incomingBasisState, coefficient) in incomingState

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
            if haskey(outgoingState, outgoingBasisState)
                outgoingState[outgoingBasisState] += opStrength * newCoefficient
            else
                outgoingState[outgoingBasisState] = opStrength * newCoefficient
            end
        end
    end
    return outgoingState
end
export ApplyOperatorChunk


"""
    ApplyOperator(operator, incomingState; tolerance)

Extends ApplyOperatorChunk() by applying a more general operator (consisting
of multiple operator chunks) on a general state.

# Examples
```jldoctest
julia> state = Dict(Bool.([1, 0]) => 1.0, Bool.([0, 1]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 0] => 1.0
  [0, 1] => -0.5

julia> operator = [("+-", [1, 2], 0.1), ("nh", [2, 1], 1.0)]
2-element Vector{Tuple{String, Vector{Int64}, Float64}}:
 ("+-", [1, 2], 0.1)
 ("nh", [2, 1], 1.0)

julia> ApplyOperator(operator, state)
Dict{BitVector, Float64} with 2 entries:
  [1, 0] => -0.05
  [0, 1] => -0.5
```
"""
function ApplyOperator(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )
    @assert !isempty(operator)
    @assert maximum([maximum(positions) for (_, positions, _) in operator]) ≤ length.(keys(incomingState))[1]

    return mergewith(+, fetch.([Threads.@spawn ApplyOperatorChunk(opType, opMembers, opStrength, copy(incomingState); tolerance=tolerance) 
                                for (opType, opMembers, opStrength) in operator])...)

    return outgoingState
end
export ApplyOperator


"""
    OperatorMatrix(basisStates, operator)

Return the matrix representation of the operator in the given basis.

# Examples
```jldoctest
julia> basis = BasisStates(2)
4-element Vector{Dict{BitVector, Float64}}:
 Dict([0, 0] => 1.0)
 Dict([0, 1] => 1.0)
 Dict([1, 0] => 1.0)
 Dict([1, 1] => 1.0)

julia> operator = [("+-", [1, 2], 0.5), ("n", [2], -1.0)]
2-element Vector{Tuple{String, Vector{Int64}, Float64}}:
 ("+-", [1, 2], 0.5)
 ("n", [2], -1.0)

julia> OperatorMatrix(basis, operator)
4×4 Matrix{Float64}:
 0.0   0.0  0.0   0.0
 0.0  -1.0  0.0   0.0
 0.0   0.5  0.0   0.0
 0.0   0.0  0.0  -1.0
```
"""
function OperatorMatrix(
        basisStates::Vector{Dict{BitVector,Float64}},
        operator::Vector{Tuple{String,Vector{Int64},Float64}};
        tolerance::Float64=1e-16,
    )
    operatorMatrix = zeros(length(basisStates), length(basisStates))
    newStates = fetch.([Threads.@spawn ApplyOperator(operator, incomingState; tolerance=tolerance)
                        for incomingState in basisStates])
    Threads.@threads for incomingIndex in findall(!isempty, newStates)
        Threads.@threads for outgoingIndex in eachindex(basisStates)
            operatorMatrix[outgoingIndex, incomingIndex] = StateOverlap(basisStates[outgoingIndex], newStates[incomingIndex])
        end
    end
    return operatorMatrix
end
export OperatorMatrix


"""
    StateOverlap(state1, state2)

Compute the inner product ⟨state1|state2⟩.

# Examples
```jldoctest
julia> state1 = Dict(Bool.([1, 0]) => 1.0, Bool.([0, 1]) => -0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 0] => 1.0
  [0, 1] => -0.5

julia> state2 = Dict(Bool.([1, 1]) => 0.5, Bool.([0, 1]) => 0.5)
Dict{BitVector, Float64} with 2 entries:
  [1, 1] => 0.5
  [0, 1] => 0.5

julia> StateOverlap(state1, state2)
-0.25
```
"""
function StateOverlap(
        state1::Dict{BitVector,Float64}, 
        state2::Dict{BitVector,Float64}
    )
    overlap = 0.
    keys2 = keys(state2)
    for (key, val) in state1
        if key ∈ keys2
            overlap += val * state2[key]
        end
    end
    return overlap
end
export StateOverlap


function ExpandIntoBasis(
        state::Dict{BitVector,Float64}, 
        basisStates::Vector{Dict{BitVector,Float64}},
    )
    coefficients = zeros(length(basisStates))
    for (index, bstate) in enumerate(basisStates)
        coefficients[index] = StateOverlap(bstate, state)
    end
    return coefficients
end
export ExpandIntoBasis


function GetSector(
        state::Dict{BitVector, Float64}, 
        symmetries::Vector{Char},
    )
    bstate = collect(keys(state))[1]
    totOcc = sum(bstate)
    magz = sum(bstate[1:2:end]) - sum(bstate[2:2:end])
    return ifelse(symmetries == ['N', 'Z'], (totOcc, magz), ifelse(symmetries == ['N'], (totOcc,), (magz,)))
end
export GetSector


"""
    roundTo(val, tolerance)

Round the provided float to the specified tolerance.

# Examples
```jldoctest
julia> fermions.roundTo(1.122323, 1e-3)
1.122
```
"""
function RoundTo(
        val::Union{Int64, Float64}, 
        tolerance::Float64
    )
    return round(val, digits=trunc(Int, -log10(tolerance)))
end
export RoundTo
