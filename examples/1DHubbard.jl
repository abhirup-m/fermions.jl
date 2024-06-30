#### Iterative diagonalisation solution of the 1D Hubbard model ####
using ProgressMeter, LinearAlgebra
include("../src/base.jl")
include("../src/correlations.jl")
include("../src/iterativeDiag.jl")

function dimerHamiltonian(U, t)
    hopping = Dict(
                   ("+-", [1, 3]) => -t,
                   ("+-", [2, 4]) => -t,
                   ("+-", [3, 1]) => -t,
                   ("+-", [4, 2]) => -t,
                  )
    correlation = Dict(
                       ("n", [1]) => -U/2,
                       ("n", [2]) => -U/2,
                       ("nn", [1, 2]) => U,
                       ("n", [3]) => -U/2,
                       ("n", [4]) => -U/2,
                       ("nn", [3, 4]) => U,
                  )
    return mergewith(+, hopping, correlation)
end


function dimerAdditionalHamiltonian(U, t, index)
    upPosition = 2 * index - 1
    hopping = Dict(
                   ("+-", [upPosition, upPosition + 2]) => -t,
                   ("+-", [upPosition + 1, upPosition + 3]) => -t,
                   ("+-", [upPosition + 2, upPosition]) => -t,
                   ("+-", [upPosition + 3, upPosition + 1]) => -t,
                   ("+-", [upPosition, upPosition - 2]) => -t,
                   ("+-", [upPosition + 1, upPosition - 1]) => -t,
                   ("+-", [upPosition - 2, upPosition]) => -t,
                   ("+-", [upPosition - 1, upPosition + 1]) => -t,
                  )
    
    correlation = Dict(
                       ("n", [upPosition]) => -U/2,
                       ("n", [upPosition + 1]) => -U/2,
                       ("nn", [upPosition, upPosition + 1]) => U
                      )
    return mergewith(+, hopping, correlation)
end

function main(numSteps, U)
    t = 1.0
    hamiltonianFamily = [dimerHamiltonian(U, t)]
    numStatesFamily = Int64[2]
    insertPosFamily = Int64[]
    insertPos = 2
    for i in 1:numSteps
        push!(hamiltonianFamily, dimerAdditionalHamiltonian(U, t, insertPos))
        push!(numStatesFamily, 2 + i)
        push!(insertPosFamily, insertPos)
        insertPos = ifelse(i % 2 ≠ 0, insertPos, insertPos + 1)
    end
    return hamiltonianFamily, numStatesFamily, insertPosFamily
end

hamiltonianFamily, numStatesFamily, insertPosFamily = main(2, 0.0)
initBasis = BasisStates(4; allowedOccupancies=0.5)
retainSize = 100
spectrumFamily = iterativeDiagonaliser(hamiltonianFamily, initBasis, numStatesFamily, insertPosFamily, retainSize)

for (numStates, pos, (E, X)) in zip(numStatesFamily, insertPosFamily, spectrumFamily[2:end])
    if numStates % 2 == 0
        continue
    end
    sector = (0.5, 0)
    groundState = X[sector][1]
    energyGs = E[sector][1]
    probe = Dict(("-", [numStates]) => 1.0)
    probeDag = Dict(("+", [numStates]) => 1.0)
    specfunc = specFunc(groundState, energyGs, E, X, probe, probeDag, (1.0, 100), 0.01)
    # function specFunc(
    #         groundState::Dict{BitVector,Float64},
    #         energyGs::Float64,
    #         eigVals::Dict{Tuple{Float64, Int64}, Vector{Float64}}, 
    #         eigVecs::Dict{Tuple{Float64, Int64}, Vector{Dict{BitVector, Float64}}},
    #         probe::Dict{Tuple{String,Vector{Int64}},Float64},
    #         probeDag::Dict{Tuple{String,Vector{Int64}},Float64},
    #         freqPoints::Tuple{Float64, Int64},
    #         broadening::Float64
    #     )
end
