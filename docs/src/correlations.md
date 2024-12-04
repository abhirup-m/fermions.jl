# Correlations - Computing Observables and Other Analyses

```@docs
GenCorrelation(
        state::Dict{BitVector,Float64},
        operator::Vector{Tuple{String,Vector{Int64},Float64}}
    )

GenCorrelation(
        state::Vector{Float64},
        operator::Matrix{Float64}
    )

ReducedDM(
        state::Dict{BitVector,Float64},
        nonTracedSites::Vector{Int64};
        reducingConfigs::Vector{BitVector}=BitVector[]
    )

PartialTraceProjectors(
        nonTracedSites::Vector{Int64};
        nonTracedConfigs::Vector{BitVector}=BitVector[]
    )

ReducedDMProjectorBased(
        state::Dict{BitVector,Float64},
        nonTracedSites::Vector{Int64};
        nonTracedConfigs::Vector{BitVector}=BitVector[]
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
    eigVecs::Vector{Vector{Float64}},
    probes::Dict{String,Matrix{Float64}},
    freqValues::Vector{Float64},
    standDev::Union{Vector{Float64}, Float64};
    degenTol::Float64=0.,
    normalise::Bool=true,
)

SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probes::Dict{String,Vector{Tuple{String,Vector{Int64},Float64}}},
    freqValues::Vector{Float64},
    basisStates::Vector{Dict{BitVector,Float64}},
    standDev::Union{Vector{Float64}, Float64};
    degenTol::Float64=0.,
    normalise::Bool=true,
)
SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probes::Dict{String,Vector{Tuple{String,Vector{Int64},Float64}}},
    freqValues::Vector{Float64},
    basisStates::Vector{Dict{BitVector,Float64}},
    standDev::Union{Vector{Float64}, Float64},
    symmetries::Vector{Char};
    degenTol::Float64=0.,
    normalise::Bool=true,
)
SpecFunc(
    eigVals::Vector{Float64},
    eigVecs::Vector{Dict{BitVector,Float64}},
    probes::Dict{String,Vector{Tuple{String,Vector{Int64},Float64}}},
    freqValues::Vector{Float64},
    basisStates::Vector{Dict{BitVector,Float64}},
    standDev::Union{Vector{Float64}, Float64},
    symmetries::Vector{Char},
    groundStateSector::Union{Tuple{Int64}, Tuple{Int64,Int64}};
    degenTol::Float64=0.,
    normalise::Bool=true,

)
SpecFunc(
    eigVals::Vector{Float64},
    eigVecMatrix::Matrix{Float64},
    probes::Dict{String,Matrix{Float64}},
    freqValues::Vector{Float64},
    standDev::Union{Vector{Float64}, Float64};
    degenTol::Float64=0.,
    normalise::Bool=true,
)
```
