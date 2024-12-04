# Base - Working with States and Operators

```@docs
BasisStates(
        numLevels::Int64, 
        totOccReq::Vector{Int64},
        magzReq::Vector{Int64},
        localCriteria::Function
)

BasisStates(
        numLevels::Int64;
        totOccReq::Union{Vector{Int64}, Int64, Nothing}=nothing,
        magzReq::Union{Vector{Int64}, Int64, Nothing}=nothing,
        localCriteria::Function=x -> true
    )

BasisStates1p(
        numLevels::Int64, 
    )

TransformBit(qubit::Bool, operator::Char)

ApplyOperatorChunk(
        opType::String,
        opMembers::Vector{Int64},
        opStrength::Float64,
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )

ApplyOperator(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )

OperatorMatrix(
        basisStates::Vector{Dict{BitVector,Float64}},
        operator::Vector{Tuple{String,Vector{Int64},Float64}}
    )

StateOverlap(
        state1::Dict{BitVector,Float64}, 
        state2::Dict{BitVector,Float64}
    )

ExpandIntoBasis(
        state::Dict{BitVector,Float64}, 
        basisStates::Vector{Dict{BitVector,Float64}},
    )

GetSector(
        state::Dict{BitVector, Float64}, 
        symmetries::Vector{Char},
    )

RoundTo(
        val::Union{Int64, Float64}, 
        tolerance::Float64
    )

PermuteSites(
        state::BitVector,
        permutation::Vector{Int64},
    )

PermuteSites(
        state::Dict{BitVector, Float64},
        permutation::Vector{Int64},
    )

Dagger(
        operator::String,
        members::Vector{Int64};
    )

Dagger(
        operator::Vector{Tuple{String,Vector{Int64},Float64}};
    )
```
