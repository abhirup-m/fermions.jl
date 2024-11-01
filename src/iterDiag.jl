using Serialization, Random, LinearAlgebra, ProgressMeter


function OrganiseOperator(
        operator::String,
        members::Vector{Int64},
    )
    if issorted(members)
        return operator, members, 1
    end

    sign = 1
    for position in length(operator):-1:1
        if sortperm(members)[position] == position
            continue
        else
            forwardSign = (-1)^length(filter(∈(('+', '-')), operator[sortperm(members)[position]+1:position]))
            backwardSign = (-1)^length(filter(∈(('+', '-')), operator[sortperm(members)[position]+1:position-1]))
            sign *= forwardSign * backwardSign
            
            bits = split(operator, "")
            tmp = bits[sortperm(members)[position]]
            bits[sortperm(members)[position]] = bits[position]
            bits[position] = tmp
            operator = join(bits)

            tmp = sort(members)[position]
            members[sortperm(members)[position]] = members[position]
            members[position] = tmp
        end
    end
    return operator, members, sign
end
export OrganiseOperator


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
        occReqSet::NTuple{2, Union{Nothing,Function}},
        magzReqSet::NTuple{2, Union{Nothing,Function}},
    )
    quantumNoReqSet = []
    for i in 1:2
        requirement = nothing
        if !isnothing(occReqSet[i]) && isnothing(magzReqSet[i])
            requirement = (q, N) -> occReqSet[i](q[1], N)
        elseif isnothing(occReqSet[i]) && !isnothing(magzReqSet[i])
            requirement = (q, N) -> magzReqSet[i](q[1], N)
        elseif !isnothing(occReqSet[i]) && !isnothing(magzReqSet[i])
            requirement = (q, N) -> occReqSet[i](q[1], N) && magzReqSet[i](q[2], N)
        end
        push!(quantumNoReqSet, requirement)
    end
    return quantumNoReqSet[1], quantumNoReqSet[2]
end


function CombineRequirements(
        occReq::Union{Nothing,Function},
        magzReq::Union{Nothing,Function},
    )
    if !isnothing(occReq) && isnothing(magzReq)
        requirement = (q, N) -> occReq(q[1], N)
    elseif isnothing(occReq) && !isnothing(magzReq)
        requirement = (q, N) -> magzReq(q[1], N)
    elseif !isnothing(occReq) && !isnothing(magzReq)
        requirement = (q, N) -> occReq(q[1], N) && magzReq(q[2], N)
    else
        requirement = (q, N) -> true
    end
    return requirement
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
    if !isnothing(quantumNos)
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
        corrQuantumNoReq::Union{Nothing,Function},
        rotation::Matrix{Float64},
        eigVals::Vector{Float64},
        maxSize::Int64,
        degenTol::Float64,
        currentSites::Vector{Int64},
        quantumNos::Union{Nothing,Vector{NTuple{1, Int64}},Vector{NTuple{2, Int64}}}, 
    )
    if isnothing(corrQuantumNoReq)
        # ensure we aren't truncating in the middle of degenerate states
        rotation = rotation[:, 1:maxSize]
        if !isnothing(quantumNos)
            quantumNos = quantumNos[1:maxSize]
        end
        eigVals = eigVals[1:maxSize]
    else
        retainBasisVectors = findall(q -> corrQuantumNoReq(q, maximum(currentSites)), quantumNos)
        otherBasisVectors = findall(q -> !corrQuantumNoReq(q, maximum(currentSites)), quantumNos)
        if length(retainBasisVectors) < maxSize
            append!(retainBasisVectors, otherBasisVectors[1:(maxSize - length(retainBasisVectors))])
        else
            retainBasisVectors = retainBasisVectors[1:maxSize]
        end

        rotation = rotation[:, retainBasisVectors]
        quantumNos = quantumNos[retainBasisVectors]
        eigVals = eigVals[retainBasisVectors]
    end
    return rotation, eigVals, quantumNos
end


function UpdateQuantumNos(
        newBasis::Vector{Dict{BitVector, Float64}},
        symmetries::Vector{Char},
        newSites::Vector{Int64}, 
        quantumNos::Union{Vector{NTuple{1, Int64}}, Vector{NTuple{2, Int64}}},
    )
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


function UpdateOldOperators(
        eigVals::Vector{Float64}, 
        identityEnv::Diagonal{Bool, Vector{Bool}}, 
        newBasket::Vector{Tuple{String, Vector{Int64}}},
        operators::Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}},
        rotation::Matrix{Float64}, 
        bondAntiSymmzer::Matrix{Float64},
        corrOperatorDict::Dict{String, Union{Nothing, Matrix{Float64}}},
        vneOperators::Vector{Tuple{String, Vector{Int64}}},
    )

    # expanded diagonal hamiltonian
    hamltMatrix = kron(diagm(eigVals), identityEnv)

    # rotate and enlarge qubit operators of current system
    for k in setdiff(keys(operators), newBasket)
        if k ∉ vneOperators
            delete!(operators, k)
        end
    end
    keySet = collect(keys(operators))

    Threads.@threads for key in keySet
        operators[key] = kron(rotation' * operators[key] * rotation, identityEnv)
    end
    bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation
    Threads.@threads for (name, corrOperator) in collect(corrOperatorDict)
        if !isnothing(corrOperator)
            corrOperatorDict[name] = kron(rotation' * corrOperator * rotation, identityEnv)
        end
    end
    return hamltMatrix, operators, bondAntiSymmzer, corrOperatorDict
end


"""Main function for iterative diagonalisation. Gives the approximate low-energy
spectrum of a hamiltonian through the following algorithm: first diagonalises a
Hamiltonian with few degrees of freedom, retains a fixed Ns number of low-energy
eigenstates, writes the Hamiltonian for a larger number of degrees of freedom in
the basis of these Ns states, then diagonalises it, etc.
"""
function IterDiag(
    hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    maxSize::Int64,
    symmetries::Vector{Char},
    correlationDefDict::Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}},
    quantumNoReq::Union{Nothing,Function},
    corrQuantumNoReq::Union{Nothing,Function},
    degenTol::Float64,
    dataDir::String,
    silent::Bool,
    specFuncNames::Vector{String},
)

    currentSites = Int64[]
    newSitesFlow = Vector{Int64}[]
    for (step, hamlt) in enumerate(hamltFlow)
        newSites = [setdiff(opMembers, currentSites) for (_, opMembers, _) in hamlt] |> V -> vcat(V...) |> unique |> sort
        push!(newSitesFlow, newSites)
        append!(currentSites, newSites)
    end

    for (i, hamlt) in enumerate(hamltFlow)
        for (j, (operatorString, members, coupling)) in enumerate(hamlt)
            newString, newMembers, sign = OrganiseOperator(operatorString, members)#, newSitesFlow[i])
            hamltFlow[i][j] = (newString, newMembers, coupling * sign)
        end
    end

    currentSites = collect(1:maximum(maximum.([opMembers for (_, opMembers, _) in hamltFlow[1]])))
    initBasis = BasisStates(maximum(currentSites))

    savePaths = SetupDataWrite(dataDir, hamltFlow, initBasis, currentSites, symmetries)

    operators, bondAntiSymmzer, hamltMatrix, newSitesFlow, create, basket = InitiateMatrices(currentSites, hamltFlow, initBasis)

    for (name, correlationDef) in correlationDefDict
        correlationMaxMember = maximum([maximum(members) for (_, members, _) in correlationDef])
        corrDefFlow = [ifelse(correlationMaxMember ∈ newSites, correlationDef, eltype(correlationDef)[]) for newSites in newSitesFlow]

        corrCreate, corrBasket = UpdateRequirements(corrDefFlow, newSitesFlow)
        create = [[c1; c2] for (c1, c2) in zip(create, corrCreate)]
        basket = [[c1; c2] for (c1, c2) in zip(basket, corrBasket)]
    end

    operators = CreateProductOperator(create[1], operators, newSitesFlow[1])

    quantumNos = InitiateQuantumNos(currentSites, symmetries, initBasis)

    corrOperatorDict = Dict{String, Union{Nothing, Matrix{Float64}}}(name => nothing for name in keys(correlationDefDict))
    if !isempty(specFuncNames)
        specFuncOperators = Dict{String, Vector{Matrix{Float64}}}(name => Matrix{Float64}[] for name in specFuncNames)
    end

    for (name, correlationDef) in correlationDefDict
        corrDefKeys = [(type, members) for (type, members, _) in correlationDef]
        if isnothing(corrOperatorDict[name]) && all(∈(keys(operators)), corrDefKeys)
            corrOperatorDict[name] = sum([coupling * operators[(type, members)] for (type, members, coupling) in correlationDef])
            if name in specFuncNames && isempty(specFuncOperators[name])
                push!(specFuncOperators[name], corrOperatorDict[name])
            end
        end
    end


    resultsDict = Dict{String, Union{Nothing, Float64}}(name => nothing for name in keys(correlationDefDict))
    @assert "energyPerSite" ∉ keys(resultsDict)
    resultsDict["energyPerSite"] = nothing

    pbar = Progress(length(hamltFlow); enabled=!silent)
    for (step, hamlt) in enumerate(hamltFlow)
        for (type, members, strength) in hamlt
            hamltMatrix += strength * operators[(type, members)]
        end
        eigVals, rotation, quantumNos = Diagonalise(hamltMatrix, quantumNos)

        if step == length(hamltFlow)
            serialize(savePaths[step], Dict("basis" => rotation,
                                            "eigVals" => eigVals,
                                            "quantumNos" => quantumNos,
                                            "currentSites" => currentSites,
                                            "newSites" => newSitesFlow[step],
                                            "bondAntiSymmzer" => bondAntiSymmzer,
                                            "results" => resultsDict,
                                           )
                     )

            resultsDict["energyPerSite"] = eigVals[1]/maximum(currentSites)

            indices = 1:length(eigVals)
            if !isnothing(corrQuantumNoReq)
                indices = findall(q -> corrQuantumNoReq(q, maximum(currentSites)), quantumNos)
            end
            finalState = rotation[:, indices[sortperm(eigVals[indices])[1]]]
            for (name, correlationDef) in correlationDefDict
                if !isnothing(corrOperatorDict[name])
                    resultsDict[name] = finalState' * corrOperatorDict[name] * finalState
                end
            end

            next!(pbar; showvalues=[("Size", size(hamltMatrix))])
            break
        end

        newBasis = BasisStates(length(newSitesFlow[step+1]))

        identityEnv = length(newSitesFlow[step+1]) == 1 ? I(2) : kron(fill(I(2), length(newSitesFlow[step+1]))...)

        saveDict = Dict("basis" => rotation,
                        "eigVals" => eigVals,
                        "quantumNos" => quantumNos,
                        "currentSites" => currentSites,
                        "newSites" => newSitesFlow[step],
                        "bondAntiSymmzer" => bondAntiSymmzer,
                        "identityEnv" => identityEnv,
                        "results" => resultsDict,
                       )

        serialize(savePaths[step], saveDict)

        if length(eigVals) > div(maxSize, 2^length(newSitesFlow[step+1]))
            rotation, eigVals, quantumNos = TruncateSpectrum(quantumNoReq, rotation, eigVals, 
                                                             div(maxSize, 2^length(newSitesFlow[step+1])), 
                                                             degenTol, currentSites, quantumNos
                                                            )
        end

        if !isnothing(quantumNos)
            quantumNos = UpdateQuantumNos(newBasis, symmetries, newSitesFlow[step+1], quantumNos)
        end

        hamltMatrix, operators, bondAntiSymmzer, corrOperator = UpdateOldOperators(eigVals, identityEnv, basket[step+1], 
                                                                                   operators, rotation, bondAntiSymmzer, 
                                                                                   corrOperatorDict, Tuple{String, Vector{Int64}}[],
                                                                                  )

        for name in filter(name -> !isempty(specFuncOperators[name]), specFuncNames)
            push!(specFuncOperators[name], corrOperatorDict[name])
        end

        # define the qbit operators for the new sites
        for site in newSitesFlow[step+1] 
            operators[("+", [site])] = kron(bondAntiSymmzer, OperatorMatrix(newBasis, [("+", [site - length(currentSites)], 1.0)]))
            operators = CreateDNH(operators, site)
        end

        operators = CreateProductOperator(create[step+1], operators, newSitesFlow[step+1])


        # construct any new correlation operators that became available at this step
        for (name, correlationDef) in correlationDefDict
            corrDefKeys = [(type, members) for (type, members, _) in correlationDef]
            if isnothing(corrOperatorDict[name]) && all(∈(keys(operators)), corrDefKeys)
                corrOperatorDict[name] = sum([coupling * operators[(type, members)] for (type, members, coupling) in correlationDef])
                if name in specFuncNames && isempty(specFuncOperators[name])
                    push!(specFuncOperators[name], corrOperatorDict[name])
                end
            end
        end

        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSitesFlow[step+1]))...)

        append!(currentSites, newSitesFlow[step+1])
        next!(pbar; showvalues=[("Size", size(hamltMatrix))])
    end
    return savePaths, resultsDict, specFuncOperators
end
export IterDiag


function IterDiag(
    hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    maxSize::Int64;
    occReq::Union{Nothing,Function}=nothing,#(o,N)->ifelse(div(N,2)-2 ≤ o ≤ div(N,2)+2, true, false), ## close to half-filling
    magzReq::Union{Nothing,Function}=nothing,#(m,N)->ifelse(m == N, true, false), ## maximally polarised states
    localCriteria::Union{Nothing,Function}=nothing,#x->x[1] + x[2] == 1, ## impurity occupancy 1
    corrOccReq::Union{Nothing,Function}=nothing,
    corrMagzReq::Union{Nothing,Function}=nothing,
    symmetries::Vector{Char}=Char[],
    degenTol::Float64=1e-10,
    dataDir::String="data-iterdiag",
    correlationDefDict::Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}=Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}(),
    vneDefDict::Dict{String, Vector{Int64}}=Dict{String, Vector{Int64}}(),
    mutInfoDefDict::Dict{String, NTuple{2,Vector{Int64}}}=Dict{String, NTuple{2,Vector{Int64}}}(),
    specFuncDefDict::Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}=Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}(),
    silent::Bool=false,
)

    exitCode = 0
    retainKeys = [k for k in keys(correlationDefDict)]
    append!(retainKeys, [k for k in keys(vneDefDict)])
    append!(retainKeys, [k for k in keys(mutInfoDefDict)])
    append!(retainKeys, [k for k in keys(specFuncDefDict)])
    push!(retainKeys, "energyPerSite")
    @assert allunique(retainKeys)

    @assert length(hamltFlow) > 1
    @assert all(∈("NS"), symmetries)
    if !isnothing(occReq)
        @assert 'N' in symmetries
    end
    if !isnothing(magzReq)
        @assert 'S' in symmetries
    end

    mutInfoToVneMap = Dict{String, Union{Nothing,NTuple{3, String}}}(name => nothing for name in keys(mutInfoDefDict))
    for (name, (sysA, sysB)) in mutInfoDefDict
        subsystemVneNames = []
        for subsystem in [sysA, sysB, vcat(sysA, sysB)]
            if subsystem ∉ values(vneDefDict)
                subsystemVneName = randstring(5)
                while subsystemVneName ∈ keys(vneDefDict)
                    subsystemVneName = randstring(5)
                end
                vneDefDict[subsystemVneName] = subsystem
            else 
                subsystemVneName = [k for (k,v) in vneDefDict if v == subsystem][1]
            end
            push!(subsystemVneNames, subsystemVneName)
        end
        mutInfoToVneMap[name] = Tuple(subsystemVneNames)
    end


    vneOperatorToCorrMap = Dict{String, Dict{Tuple{String, Vector{Int64}}, Union{Nothing,String}}}(name => Dict() for name in keys(vneDefDict))
    for (name, sites) in vneDefDict
        for sequence in Iterators.product(repeat(["+-nh"], length(sites))...)
            operator = (join(sequence), sites, 1.)
            corrName = randstring(5)
            while corrName ∈ keys(correlationDefDict)
                corrName = randstring(5)
            end
            correlationDefDict[corrName] = [operator]
            vneOperatorToCorrMap[name][(join(sequence), sites)] = corrName
        end
    end

    specFuncToCorrMap = Dict{String, String}()
    for (name, operator) in specFuncDefDict
        specFuncToCorrMap[name] = randstring(5)
        correlationDefDict[specFuncToCorrMap[name]] = operator
    end

    for (k, v) in deepcopy(correlationDefDict)
        for (j, (operatorString, members, coupling)) in enumerate(v)
            newString, newMembers, sign = OrganiseOperator(operatorString, members)#, newSitesFlow[i])
            correlationDefDict[k][j] = (newString, newMembers, coupling * sign)
        end
    end

    quantumNoReq, corrQuantumNoReq = CombineRequirements((occReq, corrOccReq), (magzReq, corrMagzReq))

    operatorMap = Dict(
                       (1,1) => "n",
                       (0,0) => "h",
                       (1,0) => "+",
                       (0,1) => "-",
                      )


    savePaths, resultsDict, specFuncOperators = IterDiag(hamltFlow, maxSize, symmetries,
                                      correlationDefDict,
                                      quantumNoReq, 
                                      corrQuantumNoReq,
                                      degenTol,
                                      dataDir,
                                      silent,
                                      collect(values(specFuncToCorrMap)),
                                     )

    if isempty(specFuncToCorrMap)
        specFuncOperators = nothing
    else
        specFuncOperators = Dict{String,Vector{Matrix{Float64}}}(name => specFuncOperators[corrName] for (name, corrName) in specFuncToCorrMap)
    end

    for (name, sites) in vneDefDict
        reducedDM = zeros(2^length(sites), 2^length(sites))
        basisVectors = collect(Iterators.product(repeat([[1, 0]], length(sites))...))
        for (i, sequence_i) in collect(enumerate(basisVectors))
            for j in eachindex(basisVectors)[i:end]
                sequence_j = basisVectors[j]
                operatorChars = prod([operatorMap[(s_i, s_j)] for (s_i, s_j) in zip(sequence_i, sequence_j)])
                corrName = vneOperatorToCorrMap[name][(operatorChars, sites)]
                matrixElement = resultsDict[corrName]
                reducedDM[i, j] = matrixElement
                reducedDM[j, i] = matrixElement'
            end
        end
        if tr(reducedDM) < 1e-10
            exitCode = 1
        else
            reducedDM /= tr(reducedDM)
        end
        rdmEigVals = filter(>(0), eigvals(Hermitian(reducedDM)))
        resultsDict[name] = -sum(rdmEigVals .* log.(rdmEigVals))
    end

    for name in keys(mutInfoDefDict)
        vneNames = mutInfoToVneMap[name]
        resultsDict[name] = resultsDict[vneNames[1]] + resultsDict[vneNames[2]] - resultsDict[vneNames[3]]
        if resultsDict[name] < -1e-10
            exitCode = 2
        end
    end

    for name in keys(resultsDict)
        if name ∉ retainKeys
            delete!(resultsDict, name)
        end
    end
    if isempty(vneDefDict) && isempty(mutInfoDefDict) && isempty(specFuncOperators)
        return savePaths, resultsDict
    elseif !(isempty(vneDefDict) && isempty(mutInfoDefDict)) && isempty(specFuncOperators)
        return savePaths, resultsDict, exitCode
    elseif isempty(vneDefDict) && isempty(mutInfoDefDict) && !isempty(specFuncOperators)
        return savePaths, resultsDict, specFuncOperators
    else
        return savePaths, resultsDict, exitCode
    end
end
export IterDiag


function IterSpecFunc(
        savePaths::Vector{String},
        specFuncOperators::Dict{String, Vector{Matrix{Float64}}},
        freqValues::Vector{Float64},
        standDev::Union{Vector{Float64}, Float64};
        occReq::Union{Nothing,Function}=nothing,
        magzReq::Union{Nothing,Function}=nothing,
        excOccReq::Union{Nothing,Function}=nothing,
        excMagzReq::Union{Nothing,Function}=nothing,
    )
    @assert issetequal(keys(specFuncOperators), ["create", "destroy"])
    quantumNoReq = CombineRequirements(occReq, magzReq)
    excQuantumNoReq = CombineRequirements(excOccReq, excMagzReq)
    totalSpecFunc = zeros(size(freqValues)...)
    
    for (i, savePath) in enumerate(savePaths[end-length(specFuncOperators["create"]):end-1])
        specFunc = zeros(size(freqValues)...)
        data = deserialize(savePath)
        basis = [collect(col) for col in eachcol(data["basis"])]
        eigVals = data["eigVals"]
        quantumNos = data["quantumNos"]
        currentSites = data["currentSites"]

        groundStateEnergy = minimum(eigVals[[quantumNoReq(q, length(currentSites)) for q in quantumNos]])
        groundStateVector = basis[sortperm(eigVals[[quantumNoReq(q, length(currentSites)) for q in quantumNos]])[1]]

        minimalStates = [groundStateVector]
        minimalEnergies = [groundStateEnergy]
        for (index, state) in enumerate(basis) 
           if !quantumNoReq(quantumNos[index], length(currentSites)) && excQuantumNoReq(quantumNos[index], length(currentSites))
               push!(minimalStates, state)
               push!(minimalEnergies, eigVals[index])
           end
       end

        specFunc = SpecFunc(minimalEnergies, minimalStates, Dict(name => specFuncOperators[name][i] for name in keys(specFuncOperators)),
                            freqValues, standDev;normalise=true)
        totalSpecFunc .+= specFunc
    end
    return totalSpecFunc
end


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
                if !(oldMembers[end] == oldMembers[1] + length(oldMembers) - 1)
                    @assert false
                end
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
    operatorMaxMember = maximum([maximum(members) for (_, members, _) in operator])
    operatorFlow = [ifelse(operatorMaxMember ∈ newSites, operator, Tuple{String,Vector{Int64},Float64}[]) for newSites in newSitesFlow]

    return UpdateRequirements(operatorFlow, newSitesFlow)
end
export UpdateRequirements
