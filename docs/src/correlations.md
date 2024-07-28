# Correlations - Computing Observables and Other Analyses

```@docs
GenCorrelation(
        state::Dict{BitVector,Float64},
        operator::Vector{Tuple{String,Vector{Int64},Float64}}
    )

ReducedDM(
        state::Dict{BitVector,Float64},
        reducingIndices::Vector{Int64};
        reducingConfigs::Vector{BitVector}=BitVector[]
    )

VonNEntropy(
        state::Dict{BitVector,Float64},
        reducingIndices::Vector{Int64};
        reducingConfigs::Vector{BitVector}=BitVector[],
        tolerance=1e-10, schmidtGap=false
    )

MutInfo(
    state::Dict{BitVector,Float64},
    reducingIndices::Tuple{Vector{Int64},Vector{Int64}};
    reducingConfigs::Tuple{Vector{BitVector},Vector{BitVector}}=(BitVector[], BitVector[])
)

TripartiteInfo(
    groundState::Dict{BitVector,Float64},
    reducingIndices::NTuple{3,Vector{Int64}};
    reducingConfigs::NTuple{3,Vector{BitVector}}=(BitVector[], BitVector[], BitVector[])
)

ThermalAverage(
    eigenStates::Vector{Dict{BitVector,Float64}},
    eigenVals::Vector{Float64},
    operator::Vector{Tuple{String,Vector{Int64},Float64}},
    invTemp::Float64,
)

SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probe::Vector{Tuple{String,Vector{Int64},Float64}},
    probeDag::Vector{Tuple{String,Vector{Int64},Float64}},
    freqArray::Vector{Float64},
    broadening::Float64;
    gsIndex::Int64=0,
)

SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probe::Vector{Tuple{String,Vector{Int64},Float64}},
    probeDag::Vector{Tuple{String,Vector{Int64},Float64}},
    freqArray::Vector{Float64},
    broadening::Float64,
    symmetries::Vector{Char};
    gsIndex::Int64=0,
)
```
