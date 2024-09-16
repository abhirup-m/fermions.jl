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


function CreateProductOperator(
        createList::Vector{Tuple{String,Vector{Int64}}},
        operators::Dict{Tuple{String,Vector{Int64}}, Matrix{Float64}},
        newSites::Vector{Int64},
    )
    uniqueList = [term for term in createList if term ∉ keys(operators)]
    if isempty(newSites)
        for (type, members) in uniqueList
            operators[(type, members)] = prod([operators[(string(o), [m])] for (o, m) in zip(type, members)])
        end
    else
        for (type, members) in uniqueList
            oldMembers = findall(∉(newSites), members)
            if !isempty(oldMembers)
                @assert oldMembers[end] == oldMembers[1] + length(oldMembers) - 1
                oldOperator = operators[(type[oldMembers], members[oldMembers])]
                leftOperator = minimum(oldMembers) > 1 ? prod([operators[(string(type[i]), [members[i]])] for i in 1:minimum(oldMembers)-1]) : 1
                rightOperator = maximum(oldMembers) < length(members) ? prod([operators[(string(type[i]), [members[i]])] for i in maximum(oldMembers)+1:length(members)]) : 1
                operators[(type, members)] = leftOperator * oldOperator * rightOperator
            else
                operators[(type, members)] = prod([operators[(string(t), [m])] for (t, m) in zip(type, members)])
            end
        end
    end
    return operators
end


function CreateDNH(
        operators::Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}},
        site::Int64,
    )
    operators[("-", [site])] = operators[("+", [site])]'
    operators[("n", [site])] = operators[("+", [site])] * operators[("-", [site])]
    operators[("h", [site])] = operators[("-", [site])] * operators[("+", [site])]
    return operators
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
    occReq::Union{Nothing,Function}=nothing,#(x,N)->ifelse(div(N,2)-2 ≤ x ≤ div(N,2)+2, true, false),
    symmetries::Vector{Char}=Char[],
    degenTol::Float64=1e-10,
    dataDir::String="data-iterdiag",
)
    @assert length(hamltFlow) > 1
    currentSites = collect(1:maximum(maximum.([opMembers for (_, opMembers, _) in hamltFlow[1]])))
    initBasis = BasisStates(maximum(currentSites))
    bondAntiSymmzer = length(currentSites) == 1 ? sigmaz : kron(fill(sigmaz, length(currentSites))...)
    create, basket = UpdateRequirements(hamltFlow)

    saveId = randstring()
    mkpath(dataDir)
    savePaths = [joinpath(dataDir, "$(saveId)-$(j)") for j in 1:length(hamltFlow)]
    push!(savePaths, joinpath(dataDir, "metadata"))
    serialize(savePaths[end], Dict("initBasis" => initBasis, "initSites" => currentSites))

    operators = Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}}()
    for site in currentSites 
        operators[("+", [site])] = OperatorMatrix(initBasis, [("+", [site], 1.0)])
        operators = CreateDNH(operators, site)
    end

    hamltMatrix = diagm(fill(0.0, 2^length(currentSites)))

    totalOcc = Int64[]
    totalNumOperator = sum([operators[("n", [i])] for i in currentSites])

    newSites = Int64[]
    operators = CreateProductOperator(create[1], operators, newSites)
    @showprogress for (step, hamlt) in enumerate(hamltFlow)
        for (type, members, strength) in hamlt
            hamltMatrix += strength * operators[(type, members)]
        end
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

        identityEnv = length(newSites) == 1 ? I(2) : kron(fill(I(2), length(newSites))...)

        serialize(savePaths[step], Dict("basis" => rotation, "eigVals" => eigVals, "newSites" => newSites, "occupancy" => totalOcc, "identityEnv" => identityEnv, "currentSites" => currentSites))

        # expanded diagonal hamiltonian
        hamltMatrix = kron(diagm(eigVals), identityEnv)

        # rotate and enlarge qubit operators of current system
        for (k,v) in collect(operators)
            if k in basket[step+1]
                operators[k] = kron(rotation' * v * rotation, identityEnv)
            else
                delete!(operators, k)
            end
        end

        # define the qbit operators for the new sites
        for site in newSites 
            operators[("+", [site])] = OperatorMatrix(newBasis, [("+", [site - length(currentSites)], 1.0)])
            operators = CreateDNH(operators, site)
        end

        newSitesOccupancy = [round(Int, col' * sum(operators[("n", [site])] for site in newSites) * col) for col in eachcol(identityEnv)]
        totalOcc = vcat([occ .+ newSitesOccupancy for occ in totalOcc]...)

        # expand new operators
        bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation
        for site in newSites 
            operators[("+", [site])] = kron(bondAntiSymmzer, operators[("+", [site])])
            operators = CreateDNH(operators, site)
        end

        operators = CreateProductOperator(create[step+1], operators, newSites)

        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSites))...)

        append!(currentSites, newSites)
    end
    return savePaths
end
export IterDiag


function IterCorrelation(savePaths::Vector{String}, corrDef::Vector{Tuple{String, Vector{Int64}, Float64}};
        occupancy::Union{Nothing, Function}=nothing)
    corrVals = Float64[]
    metadata = deserialize(savePaths[end])

    operators = Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}}()
    for site in metadata["initSites"]
        operators[("+", [site])] = OperatorMatrix(metadata["initBasis"], [("+", [site], 1.0)])
        operators = CreateDNH(operators, site)
    end
    
    corrOperator = sum([coupling .* prod([operators[(string(t), [s])] for (t,s) in zip(type, sites)]) for (type, sites, coupling) in corrDef])

    energyPerSite = []
    for (step, savePath) in enumerate(savePaths[1:end-1])
        f = deserialize(savePath)
        basis = f["basis"]
        eigVals = f["eigVals"]
        totOcc = f["occupancy"]
        currentSites = f["currentSites"]
        if isnothing(occupancy)
            gstate = basis[:,argmin(eigVals)]
            push!(energyPerSite, minimum(eigVals)/maximum(currentSites))
        else
            occupancyIndices = findall(x -> occupancy(x, maximum(currentSites)), totOcc)
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


function UpdateRequirements(hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}})
    basket = Vector{Tuple{String,Vector{Int64}}}[[] for _ in hamltFlow]
    create = Vector{Tuple{String,Vector{Int64}}}[[] for _ in hamltFlow]
    currentSites = Int64[]
    newSitesFlow = Vector{Int64}[]
    for (step, hamlt) in enumerate(hamltFlow)
        newSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamlt] |> V -> vcat(V...) |> unique |> sort
        push!(newSitesFlow, newSites)
        append!(currentSites, newSites)
    end

    for step in reverse(eachindex(hamltFlow))
        newSites = newSitesFlow[step]
        currentSites = collect(1:maximum(maximum.([opMembers for (_, opMembers, _) in hamltFlow[step]])))
        append!(create[step], [(type, members) for (type, members, _) in hamltFlow[step]])
        if step < length(hamltFlow)
            for (type, members) in basket[step+1]
                if !isempty(intersect(members, newSites))
                    push!(create[step], (type, members))
                else
                    push!(basket[step], (type, members))
                end
            end
        end
        for (type, members) in create[step]
            oldMembers = findall(∉(newSites), members)
            if !isempty(oldMembers)
                @assert oldMembers[end] == oldMembers[1] + length(oldMembers) - 1
                push!(basket[step], (type[oldMembers], members[oldMembers]))
            end
        end
        basket[step] = unique(basket[step])
    end
    return create, basket
end
