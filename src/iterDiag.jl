using Serialization, Random, LinearAlgebra, ProgressMeter


function ArrangeBasis(rotation, totalNumOperator)
    totOcc = [state' * totalNumOperator * state for state in eachcol(rotation)]
    @assert all(x -> isapprox(round(Int, x), x, atol=1e-10), totOcc)
    totOcc = round.(Int, totOcc)
    permutationOperator = zeros(length(totOcc), length(totOcc))
    for (col, row) in enumerate(sortperm(totOcc))
        permutationOperator[row, col] = 1
    end
    sectorTable = Dict{Int64, Vector{Int64}}()
    for (i, occ) in enumerate(sort(totOcc))
        if occ ∉ keys(sectorTable)
            sectorTable[occ] = [i]
        else
            push!(sectorTable[occ], i)
        end
    end
    return permutationOperator, sectorTable
end


function JoinBasis(componentBasis, eigValsTable, sectorTable, hamltMatrix)
    totalNumRows = values(componentBasis) .|> size .|> first |> sum
    totalNumCols = values(componentBasis) .|> size .|> last |> sum
    basis = zeros(totalNumRows, totalNumCols)
    eigVals = Float64[]
    rowTailPosition = 1
    colTailPosition = 1
    for k in sort(collect(keys(componentBasis)))
        sectorTable[k] = colTailPosition:colTailPosition+size(componentBasis[k])[2]-1
        basis[rowTailPosition:rowTailPosition+size(componentBasis[k])[1]-1, 
              colTailPosition:colTailPosition+size(componentBasis[k])[2]-1
             ] .= componentBasis[k]
        append!(eigVals, eigValsTable[k])
        rowTailPosition += size(componentBasis[k])[1]
        colTailPosition += size(componentBasis[k])[2]
    end
    return basis, eigVals, sectorTable
end

"""Main function for iterative diagonalisation. Gives the approximate low-energy
spectrum of a hamiltonian through the following algorithm: first diagonalises a
Hamiltonian with few degrees of freedom, retains a fixed Ns number of low-energy
eigenstates, writes the Hamiltonian for a larger number of degrees of freedom in
the basis of these Ns states, then diagonalises it, etc.
"""
function IterDiag(
    hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    maxSize::Int64;
    canonicalError::Float64=1e-10,
    symmetries::Vector{Char}=['N'],
    degenTol::Float64=1e-10,
    dataDir::String="data-iterdiag",
)
    # @assert length(hamltFlow) > 1
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
    rotation = diagm(ones(2^length(currentSites)))
    for (step, hamlt) in enumerate(hamltFlow)

        # permute basis to bring into block-diagonal form
        totalNumOperator = sum([basicMats[('n', site)] for site in currentSites])
        permutationOperator, sectorTable = ArrangeBasis(rotation, totalNumOperator)
        hamltMatrix = permutationOperator' * hamltMatrix * permutationOperator

        rotation = rotation * permutationOperator

        # rotate fermionic operators and large symmetriser
        # to be consistent with current basis
        for (k, v) in basicMats
            basicMats[k] = rotation' * v * rotation
        end
        totalNumOperator = sum([basicMats[('n', site)] for site in currentSites])
        bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation

        # get rotated+truncated basis for each symmetry sector
        componentBasis = Dict(k => zeros(length(v), length(v)) for (k, v) in sectorTable)
        eigValsTable = Dict(k => Float64[]  for k in keys(sectorTable))
        hamltMatrix += TensorProduct(hamlt, basicMats)

        Threads.@threads for (k, v) in collect(sectorTable)
            hamiltonBlock = hamltMatrix[v,v]
            F = eigen(Hermitian(hamiltonBlock))
            retainStates = ifelse(length(F.values) < maxSize, length(F.values), maxSize)
            eigValsTable[k] = real(F.values[1:retainStates])
            componentBasis[k] = real.(F.vectors[:, 1:retainStates])

            @assert all(x -> isapprox(x, k, atol=1e-10), [state' * totalNumOperator[v,v] * state
                                  for state in eachcol(F.vectors[:, 1:retainStates])]
                       )
        end

        # join small basis from each sector into complete basis
        rotation, eigVals, sectorTable = JoinBasis(componentBasis, eigValsTable, sectorTable, hamltMatrix)

        serialize(savePaths[step], Dict("operators" => basicMats, "basis" => rotation, "eigVals" => eigValsTable, "sectorTable" => sectorTable))

        if step == length(hamltFlow)
            break
        end

        additionalSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamltFlow[step+1]] |> V -> vcat(V...) |> unique |> sort
        additionalBasis = BasisStates(length(additionalSites))

        # identity matrix for the sites being added. will be `kron'ed
        # with the operators of the current sites to expand them.
        identityEnv = length(additionalSites) == 1 ? I(2) : kron(fill(I(2), length(additionalSites))...)

        # expand the basis
        rotation = kron(rotation, identityEnv)

        # expanded diagonal hamiltonian
        hamltMatrix = kron(diagm(eigVals), identityEnv)
        
        # rotate and enlarge qubit operators of current system
        for (k,v) in collect(basicMats)
            basicMats[k] = kron(v, identityEnv)
        end

        # expand the qbit operators for the additional sites into the old system
        for site in additionalSites
            basicMats[('+', site)] = kron(bondAntiSymmzer, OperatorMatrix(additionalBasis, [("+", [site - length(currentSites)], 1.0)]))
            basicMats[('-', site)] = basicMats[('+', site)]'
            basicMats[('n', site)] = basicMats[('+', site)] * basicMats[('-', site)]
            basicMats[('h', site)] = basicMats[('-', site)] * basicMats[('+', site)]
        end
        
        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(additionalSites))...)
        
        append!(currentSites, additionalSites)
    end
    return savePaths
end
export IterDiag


function IterCorrelation(savePath, corrDef, occupancy)
    f = deserialize(savePath)
    basicMats = f["operators"]
    basis = f["basis"]
    indices = f["sectorTable"][occupancy]
    groundState = basis[:,indices[1]]
    corrOperator = TensorProduct(corrDef, basicMats)
    return GenCorrelation(groundState, corrOperator)
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
