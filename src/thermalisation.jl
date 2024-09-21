function OperatorTimeEvol(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        hamiltonian::Vector{Tuple{String,Vector{Int64},Float64}},
        initState::Dict{BitVector,Float64}, 
        basisStates::Vector{Dict{BitVector,Float64}},
        deltaTime::Float64,
        numSteps::Int64;
    )
    initStateVector = ExpandIntoBasis(initState, basisStates)
    initStateVector ./= norm(initStateVector)
    operatorMatrix = convert(Matrix{ComplexF32}, OperatorMatrix(basisStates, operator))
    hamiltonianMatrix = OperatorMatrix(basisStates, hamiltonian)
    deltaUnitary = convert(Matrix{ComplexF32}, (I(size(hamiltonianMatrix)[1]) .- 1im .* hamiltonianMatrix .* deltaTime ./ 2) / (I(size(hamiltonianMatrix)[1]) .+ 1im .* hamiltonianMatrix .* deltaTime ./ 2))
    expecValueTimeEvol = zeros(numSteps)
    operatorMatrixTimeEvol = Matrix{ComplexF64}[]
    @showprogress for step in 1:numSteps
        push!(operatorMatrixTimeEvol, operatorMatrix)
        expecValueTimeEvol[step] = real(initStateVector' * operatorMatrix * initStateVector)
        operatorMatrix = deltaUnitary' * operatorMatrix * deltaUnitary
    end
    return expecValueTimeEvol, range(0, step=deltaTime, length=numSteps), operatorMatrixTimeEvol
end
export OperatorTimeEvol


function OTOC(
        operatorMatrixTimeEvol::Vector{Matrix{ComplexF64}},
        staticMatrix::Union{Matrix{ComplexF64}, Matrix{Float64}},
        hamiltonianMatrix::Matrix{Float64};
    )
    otoc = zeros(length(operatorMatrixTimeEvol))
    groundState = eigen(hamiltonianMatrix).vectors[:, 1]
    @showprogress for step in eachindex(operatorMatrixTimeEvol)
        commutator = operatorMatrixTimeEvol[step] * staticMatrix - staticMatrix * operatorMatrixTimeEvol[step]
        otoc[step] = real(groundState' 
                          * commutator' * commutator
                          * groundState
                         )
    end
    return otoc
end
export OTOC
