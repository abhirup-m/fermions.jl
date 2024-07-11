using LinearAlgebra

function Spectrum(
    hamiltonian::Vector{Tuple{String,Vector{Int64},Float64}},
    basisStates::Vector{Dict{BitVector,Float64}};
    diagElements::Vector{Float64}=Float64[],
    tolerance::Float64=1e-10,
)
    matrix = OperatorMatrix(basisStates, hamiltonian)
    keysArr = keys.(basisStates)
    valuesArr = values.(basisStates)
    if length(diagElements) != 0
        @assert length(diagElements) == length(basisStates)
        matrix += diagm(diagElements)
    end
    eigenValues, eigenStates = eigen(matrix)

    eigenVecs = [Dict{BitVector,Float64}() for _ in eigenValues]
    Threads.@threads for (j, vector) in collect(enumerate(eachrow(eigenStates)))
        for (i, c_i) in enumerate(vector)
            mergewith!(+, eigenVecs[j], Dict(keysArr[i] .=> c_i .* valuesArr[i]))
        end
        # eigenVecs[j] = mergewith(+, [Dict(keys(basisStates[i]) .=> c_i .* values(basisStates[i]))
        # for (i, c_i) in enumerate(vector)]...)
    end
    return eigenValues, eigenVecs
end
