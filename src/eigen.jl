function getSpectrum(hamiltonian::Dict{Tuple{Int64,Int64},Matrix{Float64}}; tolerance=1e-10, progressEnabled=false)
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
        eigvals[index] = sort(eigvalues)
        eigvecs[index] = [F.vectors[:, i] for i in sortperm(eigvalues)]
        if !isnothing(pbar)
            next!(pbar)
        end
    end
    if !isnothing(pbar)
        finish!(pbar)
    end
    return eigvals, eigvecs
end


function getGstate(eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}})
    gsEnergy = minimum(minimum.(values(eigvals)))
    groundStates = Vector{Float64}[]
    blocks = []
    for (block, eigvalsBlock) in eigvals
        for (i, _) in enumerate(eigvalsBlock[eigvalsBlock .≈ gsEnergy])
            push!(groundStates, eigvecs[block][i])
            push!(blocks, block)
        end
    end
    return gsEnergy, groundStates, blocks
end

