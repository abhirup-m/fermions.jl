using Serialization, Random, LinearAlgebra, ProgressMeter


function TensorProduct(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        basicMats::Dict{Tuple{Char, Int64}, Matrix{Float64}},
    )
    totalMat = zeros(size(basicMats[('+', 1)])...)
    for (opType, opMembers, opStrength) in operator
        totalMat += opStrength * prod([basicMats[(o, m)] for (o, m) in zip(opType, opMembers)]) 
    end
    return totalMat
end
export TensorProduct


"""Main function for iterative diagonalisation. Gives the approximate low-energy
spectrum of a hamiltonian through the following algorithm: first diagonalises a
Hamiltonian with few degrees of freedom, retains a fixed Ns number of low-energy
eigenstates, writes the Hamiltonian for a larger number of degrees of freedom in
the basis of these Ns states, then diagonalises it, etc.
"""
function IterDiag(
    hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    maxSize::Int64;
    occReq::Union{Nothing,Function}=nothing,#(x,N)->ifelse(div(N,2)-2 ≤ x ≤ div(N,2)+2, true, false),
    symmetries::Vector{Char}=Char[],
    degenTol::Float64=1e-10,
    dataDir::String="data-iterdiag",
)
    @assert length(hamltFlow) > 1
    currentSites = collect(1:maximum(maximum.([opMembers for (_, opMembers, _) in hamltFlow[1]])))
    initBasis = BasisStates(maximum(currentSites))
    bondAntiSymmzer = length(currentSites) == 1 ? sigmaz : kron(fill(sigmaz, length(currentSites))...)

    saveId = randstring()
    mkpath(dataDir)
    savePaths = [joinpath(dataDir, "$(saveId)-$(j)") for j in 1:length(hamltFlow)]
    push!(savePaths, joinpath(dataDir, "metadata"))
    basicMats = Dict{Tuple{Char, Int64}, Matrix{Float64}}((type, site) =>
                                                                OperatorMatrix(initBasis, [(string(type), [site], 1.0)]) 
                                                                for site in currentSites for type in ('+', '-', 'n', 'h')
                                                               )

    serialize(savePaths[end], Dict("basicMats" => basicMats, "initSites" => currentSites))

    hamltMatrix = diagm(fill(0.0, 2^length(currentSites)))

    totalOcc = Int64[]

    @showprogress for (step, hamlt) in enumerate(hamltFlow)
        hamltMatrix += TensorProduct(hamlt, basicMats)
        rotation = zeros(size(hamltMatrix)...)
        eigVals = zeros(size(hamltMatrix)[2])
        if !isempty(totalOcc) && 'N' in symmetries
            for occ in unique(totalOcc)
                indices = findall(==(occ), totalOcc)
                F = eigen(Hermitian(hamltMatrix[indices, indices]))
                rotation[indices, indices] .= F.vectors
                eigVals[indices] .= F.values
            end
            sortSequence = sortperm(eigVals)
            rotation = rotation[:, sortperm(eigVals)]
            totalOcc = totalOcc[sortperm(eigVals)]
            eigVals = sort(eigVals)
        else
            F = eigen(Hermitian(hamltMatrix))
            rotation .= F.vectors
            eigVals .= F.values
        end

        if isempty(totalOcc)
            totalNumOperator = sum([basicMats[('n', i)] for i in currentSites])
            totalOcc = [round(Int, col' * totalNumOperator * col) for col in eachcol(rotation)]
        end

        if isnothing(occReq) && length(eigVals) > maxSize
            # ensure we aren't truncating in the middle of degenerate states
            rotation = rotation[:, eigVals .≤ eigVals[maxSize] + abs(eigVals[maxSize]) * degenTol]
            totalOcc = totalOcc[eigVals .≤ eigVals[maxSize] + abs(eigVals[maxSize]) * degenTol]
            eigVals = eigVals[eigVals .≤ eigVals[maxSize] + abs(eigVals[maxSize]) * degenTol]
        end
        if !isnothing(occReq) && length(eigVals) > maxSize
            retainBasisVectors = []
            for reqOcc in unique(totalOcc)[findall(x -> occReq(x, maximum(currentSites)), unique(totalOcc))]
                indices = findall(==(reqOcc), totalOcc)
                if length(indices) > maxSize
                    indices = indices[1:maxSize]
                end
                append!(retainBasisVectors, indices)
            end
            # if length(retainBasisVectors) > maxSize
            #     retainBasisVectors = retainBasisVectors[1:maxSize]
            # end
            rotation = rotation[:, retainBasisVectors]
            totalOcc = totalOcc[retainBasisVectors]
            eigVals = eigVals[retainBasisVectors]
        end

        if step == length(hamltFlow)
            serialize(savePaths[step], Dict("basis" => rotation, "eigVals" => eigVals, "occupancy" => totalOcc, "currentSites" => currentSites))
            break
        end

        newSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamltFlow[step+1]] |> V -> vcat(V...) |> unique |> sort
        newBasis = BasisStates(length(newSites))

        # identity matrix for the sites being added. will be `kron'ed
        # with the operators of the current sites to expand them.
        identityEnv = length(newSites) == 1 ? I(2) : kron(fill(I(2), length(newSites))...)

        serialize(savePaths[step], Dict("basis" => rotation, "eigVals" => eigVals, "newSites" => newSites, "occupancy" => totalOcc, "identityEnv" => identityEnv, "currentSites" => currentSites))

        # expanded diagonal hamiltonian
        hamltMatrix = kron(diagm(eigVals), identityEnv)

        # rotate and enlarge qubit operators of current system
        for (k,v) in collect(basicMats)
            basicMats[k] = kron(rotation' * v * rotation, identityEnv)
        end

        # define the qbit operators for the new sites
        for site in newSites for type in ['+', '-', 'n', 'h']
                basicMats[(type, site)] = OperatorMatrix(newBasis, [(string(type), [site - length(currentSites)], 1.0)])
        end end

        newSitesOccupancy = [round(Int, col' * sum(basicMats[('n', site)] for site in newSites) * col) for col in eachcol(identityEnv)]
        totalOcc = vcat([occ .+ newSitesOccupancy for occ in totalOcc]...)

        # rotate the antisymmetrizer matrix
        bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation

        # expand new operators
        for site in newSites for type in ['+', '-', 'n', 'h']
                basicMats[(type, site)] = kron(bondAntiSymmzer, basicMats[(type, site)])
        end end


        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSites))...)

        append!(currentSites, newSites)
    end
    return savePaths
end
export IterDiag


function IterCorrelation(savePaths::Vector{String}, corrDef::Vector{Tuple{String, Vector{Int64}, Float64}};
        occupancy::Int64=-1)
    corrVals = Float64[]
    participatingMembers = unique(vcat([m for (_, m, _) in corrDef]...))
    metadata = deserialize(savePaths[end])
    basicMats = metadata["basicMats"]
    corrOperator = TensorProduct(corrDef, basicMats)
    energyPerSite = []
    for (step, savePath) in enumerate(savePaths[1:end-1])
        f = deserialize(savePath)
        basis = f["basis"]
        eigVals = f["eigVals"]
        totOcc = f["occupancy"]
        currentSites = f["currentSites"]
        if occupancy < 0
            gstate = basis[:,argmin(eigVals)]
            push!(energyPerSite, minimum(eigVals)/maximum(currentSites))
        else
            occupancyIndices = findall(==(occupancy), totOcc)
            gstate = basis[:, occupancyIndices][:,argmin(eigVals[occupancyIndices])]
            push!(energyPerSite, minimum(eigVals[occupancyIndices])/maximum(currentSites))
        end
        push!(corrVals, gstate' * corrOperator * gstate)
        if step < length(savePaths) - 1
            identityEnv = f["identityEnv"]
            corrOperator = kron(basis' * corrOperator * basis, identityEnv)
        end
    end
    return corrVals, energyPerSite
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
