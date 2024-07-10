using LinearAlgebra

function Spectrum(
        hamiltonian::Vector{Tuple{String, Vector{Int64}, Float64}}, 
        basisStates::Vector{Dict{BitVector, Float64}};
        diagElements::Vector{Float64}=Float64[],
        tolerance::Float64=1e-10,
    )
    matrix = OperatorMatrix(basisStates, hamiltonian)

    if length(diagElements) != 0
        @assert length(diagElements) == length(basisStates)
        matrix += diagm(diagElements)
    end
    eigenValues, eigenStates = eigen(matrix)

    eigenVecs = [Dict{BitVector, Float64}() for _ in eigenValues]
    for (j, vector) in enumerate(eachrow(eigenStates))

        for (i, c_i) in enumerate(vector)
            if c_i == 0
                continue
            end
            mergewith!(+, eigenVecs[j], Dict(k => c_i * v for (k, v) in basisStates[i]))

        end
    end
    return eigenValues, eigenVecs
end
