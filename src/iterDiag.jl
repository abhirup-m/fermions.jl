using Serialization, Random, LinearAlgebra, ProgressMeter


function OrganiseOperator(
        operator::String,
        members::Vector{Int64},
        #=newMembers::Vector{Int64},=#
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
        quantumNoReq,
        rotation,
        eigVals,
        maxSize,
        degenTol,
        currentSites,
        quantumNos, 
    )
    if isnothing(quantumNoReq)
        # ensure we aren't truncating in the middle of degenerate states
        rotation = rotation[:, 1:maxSize]
        if !isnothing(quantumNos)
            quantumNos = quantumNos[1:maxSize]
        end
        eigVals = eigVals[1:maxSize]
    elseif !isnothing(quantumNoReq)
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
    )

    # expanded diagonal hamiltonian
    hamltMatrix = kron(diagm(eigVals), identityEnv)
    # rotate and enlarge qubit operators of current system
    for k in setdiff(keys(operators), newBasket)
        delete!(operators, k)
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
    maxSize::Int64;
    occReq::Union{Nothing,Function}=nothing,#(o,N)->ifelse(div(N,2)-2 ≤ o ≤ div(N,2)+2, true, false), ## close to half-filling
    magzReq::Union{Nothing,Function}=nothing,#(m,N)->ifelse(m == N, true, false), ## maximally polarised states
    localCriteria::Union{Nothing,Function}=nothing,#x->x[1] + x[2] == 1, ## impurity occupancy 1
    corrOccReq::Union{Nothing,Function}=nothing,
    corrMagzReq::Union{Nothing,Function}=nothing,
    symmetries::Vector{Char}=Char[],
    degenTol::Float64=1e-10,
    dataDir::String="data-iterdiag",
    correlationDefDict::Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}=Dict{String, Tuple{String, Vector{Int64}, Float64}[]}(),
    silent::Bool=false,
)
    @assert length(hamltFlow) > 1
    @assert all(∈("NS"), symmetries)
    if !isnothing(occReq)
        @assert 'N' in symmetries
    end
    if !isnothing(magzReq)
        @assert 'S' in symmetries
    end

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

    for (k, v) in correlationDefDict
        for (j, (operatorString, members, coupling)) in enumerate(v)
            newString, newMembers, sign = OrganiseOperator(operatorString, members)#, newSitesFlow[i])
            correlationDefDict[k][j] = (newString, newMembers, coupling * sign)
        end
    end

    quantumNoReq, corrQuantumNoReq = CombineRequirements((occReq, corrOccReq), (magzReq, corrMagzReq))

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
    resultsDict = Dict{String, Vector{Float64}}(name => Float64[] for name in keys(correlationDefDict))
    @assert "energyPerSite" ∉ keys(resultsDict)
    resultsDict["energyPerSite"] = Float64[]

    pbar = Progress(length(hamltFlow); enabled=!silent)#, showvalues=[("Size", size(hamltMatrix))])
    for (step, hamlt) in enumerate(hamltFlow)
        for (type, members, strength) in hamlt
            hamltMatrix += strength * operators[(type, members)]
        end
        eigVals, rotation, quantumNos = Diagonalise(hamltMatrix, quantumNos)
        push!(resultsDict["energyPerSite"], eigVals[1]/maximum(currentSites))

        for (name, correlationDef) in collect(correlationDefDict)
            if isnothing(corrOperatorDict[name]) && all(∈(keys(operators)), [(type, members) for (type, members, _) in correlationDef])
                corrOperatorDict[name] = sum([coupling * operators[(type, members)] for (type, members, coupling) in correlationDef])
            end
            if !isnothing(corrOperatorDict[name])
                indices = 1:length(eigVals)
                if !isnothing(corrQuantumNoReq)
                    indices = findall(q -> corrQuantumNoReq(q, maximum(currentSites)), quantumNos)
                end
                gstate = rotation[:, indices[1]]
                push!(resultsDict[name], gstate' * corrOperatorDict[name] * gstate)
            end
        end

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
            next!(pbar)
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
                                        "results" => resultsDict,
                                       )
                 )

        if !isnothing(quantumNos)
            quantumNos = UpdateQuantumNos(newBasis, symmetries, newSitesFlow[step+1], quantumNos)
        end

        hamltMatrix, operators, bondAntiSymmzer, corrOperator = UpdateOldOperators(eigVals, identityEnv, basket[step+1], operators, rotation, bondAntiSymmzer, corrOperatorDict)

        # define the qbit operators for the new sites
        for site in newSitesFlow[step+1] 
            operators[("+", [site])] = kron(bondAntiSymmzer, OperatorMatrix(newBasis, [("+", [site - length(currentSites)], 1.0)]))
            operators = CreateDNH(operators, site)
        end

        operators = CreateProductOperator(create[step+1], operators, newSitesFlow[step+1])

        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSitesFlow[step+1]))...)

        append!(currentSites, newSitesFlow[step+1])
        next!(pbar; showvalues=[("Size", size(hamltMatrix))])
    end
    return savePaths, resultsDict
end
export IterDiag


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
