using LinearAlgebra
function getSpectrum(hamiltonian::Dict{Tuple{Int64,Int64},Matrix{Float64}}; tolerance=1e-16, progressEnabled=false)
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
        @inbounds eigvals[index] = sort(eigvalues)
        @inbounds eigvecs[index] = [F.vectors[:, i] for i in sortperm(eigvalues)]
        if !isnothing(pbar)
            next!(pbar)
        end
    end
    if !isnothing(pbar)
        finish!(pbar)
    end
    return eigvals, eigvecs
end


function getGstate(hamiltonianMatrix::Matrix{Float64}; tolerance=1e-16)
    roundedMatrix = round.(hamiltonianMatrix, digits=trunc(Int, -log10(tolerance)))
    F = eigen(Hermitian(roundedMatrix))
    eigvalues = F.values
    groundStates = Vector{Float64}[]
    for (i, E) in enumerate(eigvalues)
        if E ≈ minimum(eigvalues)
            push!(groundStates, F.vectors[:, i])
        end
    end
    return groundStates
end


function getGstate(basisStates::Dict{Tuple{Int64, Int64}, Vector{BitVector}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}})
    gsEnergy = minimum(minimum.(values(eigvals)))
    groundStates = Dict{BitVector, Float64}[]
    for (basis, eigval, eigvec) in zip(basisStates, eigvals, eigvecs)
        for vec in eigvec[eigval .≈ gsEnergy]
            push!(groundStates, Dict(zip(basis, vec)))
        end
    end
    return groundStates
end


function reducedDM(basisStates::Vector{BitArray}, groundState::Vector{Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[])
    nonReducingIndices = setdiff(1:length(basisStates[1]), reducingIndices)
    if length(reducingConfigs) == 0
        reducingConfigs = [collect(config) for config in Iterators.product(fill([1, 0], length(reducingIndices))...)]
    end
    nonReducingConfigs = vec(collect(Iterators.product(fill([1, 0], length(nonReducingIndices))...)))
    reducingCoeffs = Dict{BitVector, Dict{BitVector, Float64}}(BitVector(config) => Dict(BitVector(configPrime) => 0.0 
                                                                        for configPrime in nonReducingConfigs) for config in reducingConfigs)
    for (state, coeff) in zip(basisStates, groundState)
        if state[reducingIndices] ∈ reducingConfigs
            reducingCoeffs[state[reducingIndices]][state[nonReducingIndices]] += coeff
        end
    end
    reducedDMatrix = zeros(length(reducingConfigs), length(reducingConfigs))
    for ((i1, c1), (i2, c2)) in Iterators.product(enumerate(reducingConfigs), enumerate(reducingConfigs))
        reducedDMatrix[i1, i2] = sum(values(merge(*, reducingCoeffs[collect(c1)], reducingCoeffs[collect(c2)])))
    end
    return reducedDMatrix
end
