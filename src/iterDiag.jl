using Serialization, Random, LinearAlgebra, ProgressMeter


function TensorProduct(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        basicMats::Dict{Tuple{Char, Int64}, Matrix{Float64}},
    )
    totalMat = sum([
                    (opStrength * 
                     prod([basicMats[(o, m)] for (o, m) in zip(opType, opMembers)]) 
                    )
                   for (opType, opMembers, opStrength) in operator]
                  )
    return totalMat
end
export TensorProduct


function ArrangeBasis(rotation::Matrix{Float64}, totalOccCurrent::Vector{Int64}, totalOccAdditional::Vector{Int64})
    totalOccUpdated = vcat([occ1 .+ totalOccAdditional for occ1 in totalOccCurrent]...)
    @assert all(x -> isapprox(x, Int(x)), round.(totalOccUpdated, digits=10))
    totalOccUpdated = round.(Int, totalOccUpdated)
    permutationOperator = zeros(length(totalOccUpdated), length(totalOccUpdated))
    for (col, row) in enumerate(sortperm(totalOccUpdated))
        permutationOperator[row, col] = 1
    end
    sectorTable = Dict{Int64, Vector{Int64}}()
    for (i, totalOcc) in enumerate(sort(totalOccUpdated))
        if totalOcc ∉ keys(sectorTable)
            sectorTable[totalOcc] = [i]
        else
            push!(sectorTable[totalOcc], i)
        end
    end
    return permutationOperator, sectorTable, totalOccUpdated
end


function JoinBasis(componentBasis, eigValsTable, numSectors::Int64)
    totalNumRows = values(componentBasis) .|> size .|> first |> sum
    totalNumCols = values(componentBasis) .|> size .|> last |> sum
    basis = zeros(totalNumRows, totalNumCols)
    totalOccCurrent = Int64[]
    eigVals = Float64[]
    rowTailPosition = 1
    colTailPosition = 1
    sectorTable = Dict{Int64, Vector{Int64}}()

    minEigVals = Dict(k => minimum(v) for (k, v) in eigValsTable)
    if numSectors > 0
        minEigValKeep = sort(collect(values(minEigVals)))[minimum((length(minEigVals), numSectors))]
    else
        minEigValKeep = maximum(collect(values(minEigVals)))
    end
    retainSectors = []
    for k in sort(collect(keys(componentBasis)))
        sectorTable[k] = colTailPosition:colTailPosition+size(componentBasis[k])[2]-1
        basis[rowTailPosition:rowTailPosition+size(componentBasis[k])[1]-1, 
              colTailPosition:colTailPosition+size(componentBasis[k])[2]-1
             ] .= componentBasis[k]
        if minEigVals[k] ≤ minEigValKeep
            append!(eigVals, eigValsTable[k])
            append!(totalOccCurrent, repeat([k], length(eigValsTable[k])))
            append!(retainSectors, collect(colTailPosition:colTailPosition+size(componentBasis[k])[2]-1))
        end
        rowTailPosition += size(componentBasis[k])[1]
        colTailPosition += size(componentBasis[k])[2]
    end

    basis = basis[:, retainSectors]
    return basis, eigVals, sectorTable, totalOccCurrent
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
    push!(savePaths, joinpath(dataDir, "metadata"))
    basicMats = Dict{Tuple{Char, Int64}, Matrix{Float64}}((type, site) =>
                                                                OperatorMatrix(initBasis, [(string(type), [site], 1.0)]) 
                                                                for site in currentSites for type in ('+', '-', 'n', 'h')
                                                               )

    serialize(savePaths[end], Dict("basicMats" => basicMats, "initSites" => currentSites))

    hamltMatrix = diagm(fill(0.0, 2^length(currentSites)))

    @showprogress for (step, hamlt) in enumerate(hamltFlow)
        hamltMatrix += TensorProduct(hamlt, basicMats)
        F = eigen(Hermitian(hamltMatrix))
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


        if step == length(hamltFlow)
            serialize(savePaths[step], Dict("basis" => rotation, "eigVals" => F.values))
            break
        end

        newSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamltFlow[step+1]] |> V -> vcat(V...) |> unique |> sort
        newBasis = BasisStates(length(newSites))

        serialize(savePaths[step], Dict("basis" => rotation, "eigVals" => F.values, "newSites" => newSites))

        # identity matrix for the sites being added. will be `kron'ed
        # with the operators of the current sites to expand them.
        identityEnv = length(newSites) == 1 ? I(2) : kron(fill(I(2), length(newSites))...)

        # expanded diagonal hamiltonian
        hamltMatrix = kron(diagm(F.values[1:retainStates]), identityEnv)

        # rotate and enlarge qubit operators of current system
        for (k,v) in collect(basicMats)
            basicMats[k] = kron(rotation' * v * rotation, identityEnv)
        end
        # @time map!(x -> kron(rotation' * x * rotation, identityEnv), values(basicMats))

        # rotate the antisymmetrizer matrix
        bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation

        # absorb the qbit operators for the new sites
        for site in newSites for type in ['+', '-', 'n', 'h']
                basicMats[(type, site)] = kron(bondAntiSymmzer, OperatorMatrix(newBasis, [(string(type), [site - length(currentSites)], 1.0)]))
        end end

        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSites))...)

        append!(currentSites, newSites)
    end
    return savePaths
end
export IterDiag

function IterDiagSymm(
    hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    maxSize::Int64;
    numSectors::Int64=0,
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

    totalNumOperator = sum([basicMats[('n', site)] for site in currentSites])
    totalOccCurrent = [round(Int, state' * totalNumOperator * state) for state in eachcol(rotation)]
    totalOccAdditional = [0]
    numAdditional = 0
    additionalBasis = []
    basicMatsAdditional = []
    for (step, hamlt) in enumerate(hamltFlow)

        # permute basis to bring into block-diagonal form
        permutationOperator, sectorTable, totalOccCurrent = ArrangeBasis(rotation, totalOccCurrent, totalOccAdditional)
        hamltMatrix = permutationOperator' * hamltMatrix * permutationOperator

        rotation = rotation * permutationOperator

        # rotate fermionic operators and large symmetriser
        # to be consistent with current basis
        for (k, v) in collect(basicMats)
            basicMats[k] = rotation' * v * rotation
        end
        bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation

        # get rotated+truncated basis for each symmetry sector
        componentBasis = Dict(k => zeros(length(v), length(v)) for (k, v) in sectorTable)
        eigValsTable = Dict(k => Float64[]  for k in keys(sectorTable))
        hamltMatrix += TensorProduct(hamlt, basicMats)

        for (k, v) in collect(sectorTable)
            hamiltonBlock = hamltMatrix[v,v]
            F = eigen(Hermitian(hamiltonBlock))
            retainStates = ifelse(length(F.values) < maxSize, length(F.values), maxSize)
            eigValsTable[k] = real(F.values[1:retainStates])
            componentBasis[k] = real.(F.vectors[:, 1:retainStates])
        end

        # join small basis from each sector into complete basis
        rotation, eigVals, sectorTable, totalOccCurrent = JoinBasis(componentBasis, eigValsTable, numSectors)

        serialize(savePaths[step], Dict("operators" => basicMats, "basis" => rotation, "eigvals" => eigValsTable, "sectorTable" => sectorTable))

        if step == length(hamltFlow)
            break
        end

        additionalSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamltFlow[step+1]] |> V -> vcat(V...) |> unique |> sort
        # identity matrix for the sites being added. will be `kron'ed
        # with the operators of the current sites to expand them.
        identityEnv = length(additionalSites) == 1 ? I(2) : kron(fill(I(2), length(additionalSites))...)
        
        if length(additionalSites) != numAdditional
            additionalBasis = BasisStates(length(additionalSites))
            numAdditional = length(additionalSites)

            basicMatsAdditional = Dict{Tuple{Char, Int64}, Matrix{Float64}}((type, site) => OperatorMatrix(additionalBasis, [(string(type), [site], 1.0)]) 
                                                                            for (site, _) in enumerate(additionalSites) for type in ('+', '-', 'n', 'h'))
            totalNumAdditional = sum([basicMatsAdditional[('n', site)] for (site, _) in enumerate(additionalSites)])
            totalOccAdditional = [round(Int, state' * totalNumAdditional * state) for state in eachcol(identityEnv)]
        end

        # expand the basis
        rotation = kron(rotation, identityEnv)

        # expanded diagonal hamiltonian
        hamltMatrix = kron(diagm(eigVals), identityEnv)
        
        # rotate and enlarge qubit operators of current system
        Threads.@threads for (k,v) in collect(basicMats)
            basicMats[k] = kron(v, identityEnv)
        end

        # expand the qbit operators for the additional sites into the old system
        Threads.@threads for (i, site) in collect(enumerate(additionalSites))
            basicMats[('+', site)] = kron(bondAntiSymmzer, basicMatsAdditional[('+', i)])
            basicMats[('-', site)] = basicMats[('+', site)]'
            basicMats[('n', site)] = basicMats[('+', site)] * basicMats[('-', site)]
            basicMats[('h', site)] = basicMats[('-', site)] * basicMats[('+', site)]
        end
        
        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(additionalSites))...)
        
        append!(currentSites, additionalSites)

        totalNumOperator = sum([basicMats[('n', site)] for site in currentSites])
    end
    return savePaths
end
export IterDiagSymm


function IterCorrelation(savePaths::Vector{String}, corrDef::Vector{Tuple{String, Vector{Int64}, Float64}})
    corrVals = Float64[]
    participatingMembers = unique(vcat([m for (_, m, _) in corrDef]...))
    metadata = deserialize(savePaths[end])
    basicMats = metadata["basicMats"]
    corrOperator = TensorProduct(corrDef, basicMats)
    for (step, savePath) in enumerate(savePaths[1:end-1])
        f = deserialize(savePath)
        basis = f["basis"]
        eigVals = f["eigVals"]
        gstate = basis[:,argmin(eigVals)]
        push!(corrVals, gstate' * corrOperator * gstate)
        if step < length(savePaths) - 1
            newSites = f["newSites"]
            identityEnv = length(newSites) == 1 ? I(2) : kron(fill(I(2), length(newSites))...)
            corrOperator = kron(basis' * corrOperator * basis, identityEnv)
        end
    end
    return corrVals
end
export IterCorrelation


function IterCorrelation(savePath::String, corrDef::Vector{Tuple{String, Vector{Int64}, Float64}}, occupancy::Float64)
    f = deserialize(savePath)
    basicMats = f["operators"]
    basis = f["basis"]
    indices = f["sectorTable"][occupancy]
    eigVals = f["eigvals"][occupancy]
    groundState = basis[:,indices[argmin(eigVals)]]
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
