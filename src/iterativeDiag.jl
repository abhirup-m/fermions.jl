using IterTools

"""Expands the basis to accomodate new 1-particle states by tacking 
product states on to existing states. For eg. |10> + |01> -> (|10> + |01>)⊗|1>.
"""
function ExpandBasis(
        basis::Vector{Dict{BitVector, Float64}}, 
        sector::Tuple{Int64, Int64}, 
        eigvals::Vector{Float64}, 
        numNew::Int64,
    )
    @assert length(basis) == length(eigvals)
    expandedBasis = Dict{typeof(sector), typeof(basis)}()
    diagElements = Dict{typeof(sector), Vector{Float64}}()

    Threads.@threads for newComb in collect(product(repeat([[0, 1]], 2*numNew)...))
        newSector = sector .+ (sum(newComb), 
                               sum(newComb[1:2:end]) - sum(newComb[2:2:end]))
        expandedBasis[newSector] = [Dict(vcat.(keys(stateDict), newComb) .=> values(stateDict)) for stateDict in basis]
        diagElements[newSector] = eigvals
    end
    return expandedBasis, diagElements
end


function ClassifyBasis(basisStates)
    classifiedBasis = Dict{Tuple{Int64, Int64}, typeof(basisStates)}()

    # Multhreading doesn't help here, already tried!
    for stateDict in basisStates
        bstate = collect(keys(stateDict))[1]
        totOcc = sum(bstate)
        totMagzArr = sum(@view bstate[1:2:end]) - sum(@view bstate[2:2:end])
        if (totOcc, totMagzArr) ∉ keys(classifiedBasis)
            classifiedBasis[(totOcc, totMagzArr)] = []
        end
        push!(classifiedBasis[(totOcc, totMagzArr)], stateDict)
    end

    return classifiedBasis
end


"""Main function for iterative diagonalisation. Gives the approximate low-energy
spectrum of a hamiltonian through the following algorithm: first diagonalises a
Hamiltonian with few degrees of freedom, retains a fixed Ns number of low-energy
eigenstates, writes the Hamiltonian for a larger number of degrees of freedom in
the basis of these Ns states, then diagonalises it, etc.
"""
function iterativeDiagonaliser(
        hamiltonianFamily::Vector{Dict{Tuple{String,Vector{Int64}},Float64}},
        initBasis::Dict{Tuple{Int64, Int64}, Vector{Dict{BitVector, Float64}}},
        numStatesFamily::Vector{Int64},
        retainSize::Int64;
        occCriteria::Function = (x,y) -> true,
        magzCriteria::Function = x -> true,
        keepfrom::String="IR",
        tolerance::Float64=1e-12
    )
    @assert length(hamiltonianFamily) == length(numStatesFamily)

    # stores flow of the spectrum
    spectrumFamily = []

    # stores the flowing basis states
    basisStates = initBasis

    # stores the current eigenvalues
    diagElements = Dict{Tuple{Int64, Int64}, Vector{Float64}}(k => zeros(length(v)) for (k, v) in basisStates)

    # loop over the Hamiltonians H_N.
    @showprogress desc="iterative diag." for (i, hamiltonian) in enumerate(hamiltonianFamily)

        # at Nth step on bringing sites 'q' into the Hamiltonian,
        # full hamiltonian H_N is H_N = H_{N-1} + H_{N-1,q} + H_q,
        # where H_q is the terms acting purely on q, while H_{N-1,q}
        # is the set of terms that couple 'q' with the N-1 states.

        println("$i. Creating matrix.")
        # matrix for the H_{N-1,q} + H_q part.
        hamiltonianMatrix = operatorMatrix(basisStates, hamiltonian)
        # Adds the matrix elements for H_{N-1}. Since the present basis
        # is the eigenbasis of H_{N-1}, these matrix elements are the diagonal ones.
        for k in keys(hamiltonianMatrix)
            hamiltonianMatrix[k] .+= diagm(diagElements[k])
        end

        println("$i. Diagonalising matrix.")
        # diagonalise this matrix, keeping only 'retainSize' 
        # number of the new eigenstates.
        spectrum = getSpectrum(basisStates, hamiltonianMatrix; maxNum=retainSize, keepfrom=keepfrom, tolerance=tolerance)

        push!(spectrumFamily, spectrum)
        @assert keys(basisStates) == keys(spectrum[1]) == keys(spectrum[2])
        # create basis out of the eigenstates
        basisStates = spectrum[2]

        println("$i. Expanding basis.")
        if i < length(numStatesFamily)
            # expand the new basis to accomodate the new sites for the next step.
            basisStates, diagElements = expandBasis(basisStates, numStatesFamily[i+1] - numStatesFamily[i], spectrum[1]; 
                                                          occCriteria=occCriteria, magzCriteria=magzCriteria)
        end
    end
    return spectrumFamily
end
