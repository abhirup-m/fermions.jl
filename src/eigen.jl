using LinearAlgebra
function getSpectrum(hamiltonian::Dict{Tuple{Int64,Int64},Matrix{Float64}}; tolerance=1e-16, maxNum=0, progressEnabled=false)
    eigvals = Dict{Tuple{Int64,Int64},Vector{Float64}}(index => [] for index in keys(hamiltonian))
    eigvecs = Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}(index => [] for index in keys(hamiltonian))
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
        @inbounds eigvecs[index] = [round.(F.vectors[:, i], digits=trunc(Int, -log10(tolerance))) 
                                    for i in sortperm(eigvalues)[1:maxNum]]
        if !isnothing(pbar)
            next!(pbar)
        end
    end
    if !isnothing(pbar)
        finish!(pbar)
    end
    return eigvals, eigvecs
end


function getGstate(hamiltonianMatrix::Matrix{Float64}; tolerance=1e-10)
    roundedMatrix = round.(hamiltonianMatrix, digits=trunc(Int, -log10(tolerance)))
    F = eigen(Hermitian(roundedMatrix))
    eigvalues = F.values
    groundStates = Vector{Float64}[]
    for (i, E) in enumerate(eigvalues)
        if E ≈ minimum(eigvalues) atol=tolerance
            push!(groundStates, F.vectors[:, i])
        end
    end
    return groundStates
end


function getGstate(basisStates::Dict{Tuple{Int64, Int64}, Vector{BitVector}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}})
    gsEnergy = minimum(minimum.(values(eigvals)))
    groundStates = Dict{BitVector, Float64}[]
    for key in keys(basisStates)
        basis = basisStates[key]
        eigval = eigvals[key]
        eigvec = eigvecs[key]
        for vec in eigvec[eigval .≈ gsEnergy]
            @assert sum(vec .^ 2) ≈ 1
            @assert length(basis) == length(vec)
            push!(groundStates, Dict(zip(basis, vec)))
        end
    end
    return groundStates
end


function reducedDM(groundState::Dict{BitVector, Float64}, reducingIndices::Vector{Int64}; reducingConfigs::Vector{BitVector}=BitVector[])
    nonReducingIndices = setdiff(1:length(collect(keys(groundState))[1]), reducingIndices)
    if length(reducingConfigs) == 0
        reducingConfigs = [collect(config) for config in Iterators.product(fill([1, 0], length(reducingIndices))...)]
    end
    nonReducingConfigs = vec(collect(Iterators.product(fill([1, 0], length(nonReducingIndices))...)))
    reducingCoeffs = Dict{BitVector, Dict{BitVector, Float64}}(BitVector(config) => Dict(BitVector(configPrime) => 0.0 
                                                                        for configPrime in nonReducingConfigs) for config in reducingConfigs)
    for (state, coeff) in groundState
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
