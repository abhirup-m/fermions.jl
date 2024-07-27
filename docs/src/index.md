# Fermions.jl documentation

## Base.jl - Basic Methods
```@docs
BasisStates(
        numLevels::Int64, 
        totOccReq::Vector{Int64},
        magzReq::Vector{Int64},
        localCriteria::Function
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
```
