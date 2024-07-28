# Eigen - Eigenfunction Solver and Related Methods

```@docs
ClassifyBasis(
        basisStates::Vector{Dict{BitVector,Float64}},
        symmetries::Vector{Char};
        energies::Vector{Float64}=Float64[],
    )

TransformState(
        vector::Vector{Float64},
        basisStates::Vector{Dict{BitVector,Float64}};
        tolerance::Float64=1e-16,
    )

Spectrum(
    operator::Vector{Tuple{String,Vector{Int64},Float64}},
    basisStates::Vector{Dict{BitVector,Float64}};
    diagElements::Vector{Float64}=Float64[],
    tolerance::Float64=1e-16,
    )

Spectrum(
    operator::Vector{Tuple{String,Vector{Int64},Float64}},
    basisStates::Vector{Dict{BitVector,Float64}},
    symmetries::Vector{Char};
    diagElements::Vector{Float64}=Float64[],
    tolerance::Float64=1e-16,
    classify::Bool=false,
)
```
