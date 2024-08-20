using Serialization, Random, LinearAlgebra, ProgressMeter

"""Main function for iterative diagonalisation. Gives the approximate low-energy
spectrum of a hamiltonian through the following algorithm: first diagonalises a
Hamiltonian with few degrees of freedom, retains a fixed Ns number of low-energy
eigenstates, writes the Hamiltonian for a larger number of degrees of freedom in
the basis of these Ns states, then diagonalises it, etc.
"""
function IterDiag(
    hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    maxSize::Int64;
    degenTol::Float64=1e-10,
    dataDir::String="data-iterdiag",
)
    @assert length(hamltFlow) > 1
    currentSites = collect(1:maximum(maximum.([opMembers for (_, opMembers, _) in hamltFlow[1]])))
    initBasis = BasisStates(maximum(currentSites))
    bondAntiSymmzer = length(currentSites) == 1 ? sigmaz : kron(fill(sigmaz, length(currentSites))...)

    saveId = randstring()
    rm(dataDir; recursive=true, force=true)
    mkpath(dataDir)
    savePaths = [joinpath(dataDir, "$(saveId)-$(j)") for j in 1:length(hamltFlow)]
    basicMats = Dict{Tuple{Char, Int64}, Matrix{Float64}}((type, site) => OperatorMatrix(initBasis, [(string(type), [site], 1.0)]) 
                                                          for site in currentSites for type in ('+', '-', 'n', 'h'))
    hamltMatrix = diagm(fill(0.0, length(initBasis)))

    @showprogress for (step, hamlt) in enumerate(hamltFlow)
        hamltMatrix += TensorProduct(hamlt, basicMats)
        F = eigen(0.5 * (hamltMatrix + hamltMatrix'))
        retainStates = ifelse(length(F.values) < maxSize, length(F.values), maxSize)

        # ensure we aren't truncating in the middle of degenerate states
        for energy in F.values[retainStates+1:end]
            if abs(1 - F.values[retainStates] / energy) > degenTol
                break
            else
                retainStates += 1
            end
        end

        rotation = F.vectors[:, 1:retainStates]
        eigVals = F.values[1:retainStates]
        serialize(savePaths[step], Dict("operators" => basicMats, "rotation" => rotation, "eigVals" => eigVals))

        if step == length(hamltFlow)
            break
        end

        newSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamltFlow[step+1]] |> V -> vcat(V...) |> unique |> sort
        newBasis = BasisStates(length(newSites))

        # identity matrix for the sites being added. will be `kron'ed
        # with the operators of the current sites to expand them.
        identityEnv = length(newSites) == 1 ? I(2) : kron(fill(I(2), length(newSites))...)

        # expanded diagonal hamiltonian
        hamltMatrix = kron(diagm(eigVals), identityEnv)

        # rotate and enlarge qubit operators of current system
        for (k,v) in collect(basicMats)
            basicMats[k] = kron(rotation' * v * rotation, identityEnv)
        end
        
        # rotate the antisymmetrizer matrix
        bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation

        # absorb the qbit operators for the new sites
        for site in newSites for type in ('+', '-', 'n', 'h')
            basicMats[(type, site)] = kron(bondAntiSymmzer, OperatorMatrix(newBasis, [(string(type), [site - length(currentSites)], 1.0)]))
            end
        end

        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSites))...)
        
        append!(currentSites, newSites)
    end
    return savePaths
end
export IterDiag


function IterCorrelation(savePath, corrDef)
    f = deserialize(savePath)
    basicMats = f["operators"]
    rotation = f["rotation"]
    eigVals = f["eigVals"]
    eigVecs = [collect(vec) for vec in eachcol(rotation)]
    corrOperator = TensorProduct(corrDef, basicMats)
    return GenCorrelation(eigVecs[1], corrOperator)
end
export IterCorrelation


function IterSpecFunc(savePaths, probe, probeDag, freqArray, broadening)
    specfuncFlow = [0 .* freqArray for _ in savePaths]

    for (step, savePath) in enumerate(savePaths)
        f = deserialize(savePath)
        basicMats = f["operators"]
        rotation = f["rotation"]
        eigVals = f["eigVals"]
        eigVecs = [collect(vec) for vec in eachcol(rotation)]
        probeMatrix = TensorProduct(probe, basicMats)
        probeDagMatrix = TensorProduct(probeDag, basicMats)
        specfuncFlow[step] .+= SpecFunc(eigVals, rotation, probeMatrix, probeDagMatrix, freqArray, broadening)
        for j in step+1:length(specfuncFlow)
            specfuncFlow[j] = copy(specfuncFlow[step])
        end
        run(`rm $(savePath)`)
    end
    specfuncFlow[end] ./= sum(specfuncFlow[end]) * abs(freqArray[2] - freqArray[1])
    return specfuncFlow, freqArray
end
export IterSpecFunc


function TensorProduct(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        basicMats::Dict{Tuple{Char, Int64}, Matrix{Float64}},
    )
    totalMat = 0 .* collect(values(basicMats))[1]
    for (opType, opMembers, opStrength) in operator
        totalMat += opStrength * prod([basicMats[pair] for pair in zip(opType, opMembers)])
    end
    return totalMat
end
export TensorProduct
