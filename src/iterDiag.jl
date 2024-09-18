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


function CombineRequirements(
        occReq::Union{Nothing,Function},
        magzReq::Union{Nothing,Function}
    )
    quantumNoReq = nothing
    if !isnothing(occReq) && isnothing(magzReq)
        quantumNoReq = (q, N) -> occReq(q[1], N)
    elseif isnothing(occReq) && !isnothing(magzReq)
        quantumNoReq = (q, N) -> magzReq(q[1], N)
    elseif !isnothing(occReq) && !isnothing(magzReq)
        quantumNoReq = (q, N) -> occReq(q[1], N) && magzReq(q[2], N)
    end
    return quantumNoReq
end


function InitiateQuantumNos(
        currentSites::Vector{Int64}, 
        symmetries::Vector{Char},
        initBasis::Vector{Dict{BitVector, Float64}},
    )
    symmetryOperators = Vector{Tuple{String,Vector{Int64},Float64}}[]
    if 'N' in symmetries
        push!(symmetryOperators, [("n", [i], 1.) for i in currentSites])
    end
    if 'S' in symmetries
        push!(symmetryOperators, [("n", [i], (-1)^(i+1)) for i in currentSites])
    end

    if !isempty(symmetries)
        return [Tuple(round(Int, GenCorrelation(state, operator)) for operator in symmetryOperators) 
                      for state in initBasis]
    else
        return nothing
    end
end


function InitiateMatrices(currentSites, hamltFlow, initBasis)
    bondAntiSymmzer = length(currentSites) == 1 ? sigmaz : kron(fill(sigmaz, length(currentSites))...)
    create, basket, newSitesFlow = UpdateRequirements(hamltFlow)
    hamltMatrix = diagm(fill(0.0, 2^length(currentSites)))
    operators = Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}}()
    for site in currentSites 
        operators[("+", [site])] = OperatorMatrix(initBasis, [("+", [site], 1.0)])
        operators = CreateDNH(operators, site)
    end
    operators = CreateProductOperator(create[1], operators, newSitesFlow[1])
    return operators, bondAntiSymmzer, hamltMatrix, newSitesFlow, create, basket
end


function SetupDataWrite(dataDir, hamltFlow, initBasis, currentSites, symmetries)
    saveId = randstring()
    mkpath(dataDir)
    savePaths = [joinpath(dataDir, "$(saveId)-$(j)") for j in 1:length(hamltFlow)]
    push!(savePaths, joinpath(dataDir, "metadata"))
    serialize(savePaths[end], Dict("initBasis" => initBasis,
                                   "initSites" => currentSites,
                                   "symmetries" => symmetries)
             )
    return savePaths
end


function Diagonalise(
        hamltMatrix::Matrix{Float64},
        quantumNos::Union{Nothing, Vector{NTuple{1, Int64}}, Vector{NTuple{2, Int64}}},
    )
    rotation = zeros(size(hamltMatrix)...)
    eigVals = zeros(size(hamltMatrix)[2])
    @time if !isnothing(quantumNos)
        numCaptured = 0
        for quantumNo in unique(quantumNos)
            indices = findall(==(quantumNo), quantumNos)
            numCaptured += length(indices)
            F = eigen(Hermitian(hamltMatrix[indices, indices]))
            rotation[indices, indices] .= F.vectors
            eigVals[indices] .= F.values
        end
        @assert numCaptured == length(quantumNos)
        sortSequence = sortperm(eigVals)
        rotation = rotation[:, sortperm(eigVals)]
        quantumNos = quantumNos[sortperm(eigVals)]
        eigVals = sort(eigVals)
    else
        F = eigen(Hermitian(hamltMatrix))
        rotation .= F.vectors
        eigVals .= F.values
    end
    return eigVals, rotation, quantumNos
end


function TruncateSpectrum(
        quantumNoReq,
        rotation,
        eigVals,
        maxSize,
        degenTol,
        currentSites,
        quantumNos, 
    )
    if isnothing(quantumNoReq) && length(eigVals) > maxSize
        # ensure we aren't truncating in the middle of degenerate states
        rotation = rotation[:, eigVals .≤ eigVals[maxSize] + abs(eigVals[maxSize]) * degenTol]
        if !isnothing(quantumNos)
            quantumNos = quantumNos[eigVals .≤ eigVals[maxSize] + abs(eigVals[maxSize]) * degenTol]
        end
        eigVals = eigVals[eigVals .≤ eigVals[maxSize] + abs(eigVals[maxSize]) * degenTol]
    end
    if !isnothing(quantumNoReq) && length(eigVals) > maxSize
        retainBasisVectors = []
        for sector in unique(quantumNos)[findall(q -> quantumNoReq(q, maximum(currentSites)), unique(quantumNos))]
            indices = findall(==(sector), quantumNos)
            if length(indices) > maxSize
                indices = indices[1:maxSize]
            end
            append!(retainBasisVectors, indices)
        end

        rotation = rotation[:, retainBasisVectors]
        quantumNos = quantumNos[retainBasisVectors]
        eigVals = eigVals[retainBasisVectors]
    end
    return rotation, eigVals, quantumNos
end


function UpdateQuantumNos(newBasis, symmetries, newSites, quantumNos)
    newSymmetryOperators = []
    if 'N' in symmetries
        push!(newSymmetryOperators, [("n", [i], 1.) for i in eachindex(newSites)])
    end
    if 'S' in symmetries
        push!(newSymmetryOperators, [("n", [i], (-1.)^(i+1)) for i in eachindex(newSites)])
    end

    if !isempty(symmetries)
        newSitesQuantumNos = [Tuple(round(Int, GenCorrelation(state, operator)) for operator in newSymmetryOperators) 
                              for state in newBasis]
        quantumNos = vcat([[quantumNo .+ newQuantumNo for newQuantumNo in newSitesQuantumNos]
                           for quantumNo in quantumNos]...)
    end
    return quantumNos
end


function UpdateOldOperators(eigVals, identityEnv, newBasket, operators, rotation, bondAntiSymmzer)
    # expanded diagonal hamiltonian
    hamltMatrix = kron(diagm(eigVals), identityEnv)

    # rotate and enlarge qubit operators of current system
    for k in setdiff(keys(operators), newBasket)
        delete!(operators, k)
    end
    @time Threads.@threads for (k, v) in collect(operators)
        operators[k] = kron(rotation' * v * rotation, identityEnv)
    end

    bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation
    return hamltMatrix, operators, bondAntiSymmzer
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
    occReq::Union{Nothing,Function}=nothing,#(o,N)->ifelse(div(N,2)-2 ≤ o ≤ div(N,2)+2, true, false), ## close to half-filling
    magzReq::Union{Nothing,Function}=nothing,#(m,N)->ifelse(m == N, true, false), ## maximally polarised states
    symmetries::Vector{Char}=Char[],
    degenTol::Float64=1e-10,
    dataDir::String="data-iterdiag",
)
    @assert length(hamltFlow) > 1
    @assert all(∈("NS"), symmetries)
    if !isnothing(occReq)
        @assert 'N' in symmetries
    end
    if !isnothing(magzReq)
        @assert 'S' in symmetries
    end

    quantumNoReq = CombineRequirements(occReq, magzReq)

    currentSites = collect(1:maximum(maximum.([opMembers for (_, opMembers, _) in hamltFlow[1]])))
    initBasis = BasisStates(maximum(currentSites))

    savePaths = SetupDataWrite(dataDir, hamltFlow, initBasis, currentSites, symmetries)

    operators, bondAntiSymmzer, hamltMatrix, newSitesFlow, create, basket = InitiateMatrices(currentSites, hamltFlow, initBasis)
    quantumNos = InitiateQuantumNos(currentSites, symmetries, initBasis)

    energyPerSite = Float64[]
    @showprogress for (step, hamlt) in enumerate(hamltFlow)
        @time for (type, members, strength) in hamlt
            hamltMatrix += strength * operators[(type, members)]
        end
        println("Hamiltonian size = ", size(hamltMatrix))
        eigVals, rotation, quantumNos = Diagonalise(hamltMatrix, quantumNos)
        push!(energyPerSite, eigVals[1]/maximum(currentSites))

        if step == length(hamltFlow)
            serialize(savePaths[step], Dict("basis" => rotation,
                                            "eigVals" => eigVals,
                                            "quantumNos" => quantumNos,
                                            "currentSites" => currentSites,
                                            "newSites" => newSitesFlow[step],
                                            "bondAntiSymmzer" => bondAntiSymmzer,
                                           )
                     )
            break
        end

        if length(eigVals) > maxSize
            rotation, eigVals, quantumNos = TruncateSpectrum(quantumNoReq, rotation, eigVals, maxSize, degenTol, currentSites, quantumNos)
        end

        newBasis = BasisStates(length(newSitesFlow[step+1]))

        identityEnv = length(newSitesFlow[step+1]) == 1 ? I(2) : kron(fill(I(2), length(newSitesFlow[step+1]))...)

        serialize(savePaths[step], Dict("basis" => rotation,
                                        "eigVals" => eigVals,
                                        "quantumNos" => quantumNos,
                                        "currentSites" => currentSites,
                                        "newSites" => newSitesFlow[step],
                                        "bondAntiSymmzer" => bondAntiSymmzer,
                                        "identityEnv" => identityEnv,
                                       )
                 )


        quantumNos = UpdateQuantumNos(newBasis, symmetries, newSitesFlow[step+1], quantumNos)

        @time hamltMatrix, operators, bondAntiSymmzer = UpdateOldOperators(eigVals, identityEnv, basket[step+1], operators, rotation, bondAntiSymmzer)

        # define the qbit operators for the new sites
        @time Threads.@threads for site in newSitesFlow[step+1] 
            operators[("+", [site])] = kron(bondAntiSymmzer, OperatorMatrix(newBasis, [("+", [site - length(currentSites)], 1.0)]))
            operators = CreateDNH(operators, site)
        end

        @time operators = CreateProductOperator(create[step+1], operators, newSitesFlow[step+1])

        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSitesFlow[step+1]))...)

        append!(currentSites, newSitesFlow[step+1])
    end
    return savePaths, energyPerSite
end
export IterDiag


function IterCorrelation(
        savePaths::Vector{String}, 
        corrDef::Vector{Tuple{String, Vector{Int64}, Float64}};
        occReq::Union{Nothing, Function}=nothing,
        magzReq::Union{Nothing, Function}=nothing,
    )

    newSitesFlow = [deserialize(savePath)["newSites"] for savePath in savePaths[1:end-1]]
    create, basket, newSitesFlow = UpdateRequirements(corrDef, newSitesFlow)

    quantumNoReq = nothing
    if !isnothing(occReq) && isnothing(magzReq)
        quantumNoReq = (q, N) -> occReq(q[1], N)
    elseif isnothing(occReq) && !isnothing(magzReq)
        quantumNoReq = (q, N) -> magzReq(q[1], N)
    elseif !isnothing(occReq) && !isnothing(magzReq)
        quantumNoReq = (q, N) -> occReq(q[1], N) && magzReq(q[2], N)
    end
    corrVals = Float64[]
    metadata = deserialize(savePaths[end])
    symmetries = metadata["symmetries"]
    if !isnothing(occReq)
        @assert 'N' in symmetries
    end
    if !isnothing(magzReq)
        @assert 'S' in symmetries
    end

    operators = Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}}()
    for site in metadata["initSites"]
        operators[("+", [site])] = OperatorMatrix(metadata["initBasis"], [("+", [site], 1.0)])
        operators = CreateDNH(operators, site)
    end
    operators = CreateProductOperator(create[1], operators, newSitesFlow[1])

    corrOperator = nothing
    if all(∈(keys(operators)), [(type, members) for (type, members, _) in corrDef])
        corrOperator = sum([coupling * operators[(type, members)] for (type, members, coupling) in corrDef])
    end
    
    for (step, savePath) in enumerate(savePaths[1:end-1])
        f = deserialize(savePath)
        basis = f["basis"]
        eigVals = f["eigVals"]
        quantumNos = f["quantumNos"]
        currentSites = f["currentSites"]
        if !isnothing(corrOperator)
            if isnothing(quantumNoReq)
                gstate = basis[:,argmin(eigVals)]
            else
                indices = findall(q -> quantumNoReq(q, maximum(currentSites)), quantumNos)
                gstate = basis[:, indices][:,argmin(eigVals[indices])]
            end
            push!(corrVals, gstate' * corrOperator * gstate)
            if step < length(savePaths) - 1
                corrOperator = kron(basis' * corrOperator * basis, f["identityEnv"])
            end
        else
            for (k, v) in collect(operators)
                operators[k] = kron(basis' * v * basis, f["identityEnv"])
            end

            newBasis = BasisStates(length(newSitesFlow[step+1]))

            # define the qbit operators for the new sites
            for site in newSitesFlow[step+1] 
                operators[("+", [site])] = kron(basis' * f["bondAntiSymmzer"] * basis, OperatorMatrix(newBasis, [("+", [site - length(currentSites)], 1.0)]))
                operators = CreateDNH(operators, site)
            end

            operators = CreateProductOperator(create[step+1], operators, newSitesFlow[step+1])
            if all(∈(keys(operators)), [(type, members) for (type, members, _) in corrDef])
                corrOperator = sum([coupling * operators[(type, members)] for (type, members, coupling) in corrDef])
            end

        end
    end
    return corrVals
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


function UpdateRequirements(
        hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
        newSitesFlow::Vector{Vector{Int64}},
    )
    basket = Vector{Tuple{String,Vector{Int64}}}[[] for _ in hamltFlow]
    create = Vector{Tuple{String,Vector{Int64}}}[[] for _ in hamltFlow]

    for step in reverse(eachindex(hamltFlow))
        newSites = newSitesFlow[step]
        if !isempty(hamltFlow[step])
            append!(create[step], [(type, members) for (type, members, _) in hamltFlow[step]])
        end
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
    return create, basket, newSitesFlow
end
export UpdateRequirements


function UpdateRequirements(
        hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    )
    currentSites = Int64[]
    newSitesFlow = Vector{Int64}[]
    for (step, hamlt) in enumerate(hamltFlow)
        newSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamlt] |> V -> vcat(V...) |> unique |> sort
        push!(newSitesFlow, newSites)
        append!(currentSites, newSites)
    end

    return UpdateRequirements(hamltFlow, newSitesFlow)
end
export UpdateRequirements


function UpdateRequirements(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        newSitesFlow::Vector{Vector{Int64}},
    )
    currentSites = Int64[]
    operatorMaxMember = maximum([maximum(members) for (_, members, _) in operator])
    operatorFlow = [ifelse(operatorMaxMember ∈ newSites, operator, Tuple{String,Vector{Int64},Float64}[]) for newSites in newSitesFlow]

    return UpdateRequirements(operatorFlow, newSitesFlow)
end
export UpdateRequirements
