using Serialization, Random, LinearAlgebra, ProgressMeter


"""
Rearrange operator such that indices appear in ascending order.
c^†_1 c_3 c_2 c^†_4 ⟶ -c^†_1 c_2 c_3 c^†_4.
This allows us to split an operator into old parts and new parts.
For eg., if indices 1 and 2 enter at the first step and 3 and 4 
at the second step, then the sequence (1,3,2,4) cannot be split
but (1,2,3,4) can be.
"""
function OrganiseOperator(
        operator::String,
        members::Vector{Int64},
    )
    # if already sorted, return immediately
    if issorted(members)
        return operator, members, 1
    end

    sign = 1

    # calculate the sign changes associated with each index exchange.
    # The idea is to start from the end of the operator, consider each c-operator
    # and slot it into its appropriate sorted position. This incurs sign changes,
    # which we take care of by multiplying with the operator.
    for position in length(operator):-1:1
        if sortperm(members)[position] == position
            continue
        else

            sortPos = sortperm(members)[position]

            # sign change incurred when the operator at the end is teleported to its
            # correct sorted position `p'.
            if operator[position] ∈ ('+', '-')
                forwardSign = (-1)^count(∈(('+', '-')), operator[sortPos:position-1])
            else
                forwardSign = 1
            end

            # sign change incurred when the operator currently sitting at `p' is 
            # teleported to the position left vacant by the previous teleportation.
            if operator[sortPos] ∈ ('+', '-')
                backwardSign = (-1)^count(∈(('+', '-')), operator[sortPos+1:position-1])
            else
                backwardSign = 1
            end

            sign *= forwardSign * backwardSign
            
            # perform the actual teleportation. This requires converting the string
            # into a vector, because strings are immutable.
            bits = split(operator, "")

            tmp = bits[sortperm(members)[position]]
            bits[sortperm(members)[position]] = bits[position]
            bits[position] = tmp

            tmp = sort(members)[position]
            members[sortperm(members)[position]] = members[position]
            members[position] = tmp

            # convert the teleported array back into a string
            operator = join(bits)
        end
    end
    return operator, members, sign
end
export OrganiseOperator


"""
Create the tensor product O_i O_j O_k ... O_m, where the smaller operators
{O_i} are already given as matrix representations. If all the indices are
already present in the system, then it is a simple tensor product. If one
or more are new indices, then first track the composite operators for the
old and new sets, then tensor multiply these two composite operators.
"""
function CreateProductOperator(
        productOperatorDef::Tuple{String,Vector{Int64}},
        operators::Dict{Tuple{String,Vector{Int64}}, Matrix{Float64}};
        newSites::Vector{Int64}=Int64[],
    )

    # if the requested operator is already part of the operator set,
    # then nothing to do.
    if productOperatorDef ∈ keys(operators)
        return operators
    end

    # type is the string of c and cdags. members is the list of indices
    # that specify where these c and cdags will act.
    type, members = productOperatorDef
    @assert issorted(members)

    # if there are no new sites, then it is a simple tensor product.
    if isempty(newSites)
        operators[(type, members)] = prod([operators[(string(o), [m])] for (o, m) in zip(type, members)])
    else

        # if there are new sites, first extract the positions of the
        # old and new sites.
        oldIndices = findall(∉(newSites), members)
        newIndices = findall(∈(newSites), members)

        # if it isn't all new indices, caculate the old composite
        # operator.
        oldOperator = 1
        if !isempty(oldIndices)
            @assert oldIndices == 1:length(oldIndices)
            oldOperator = operators[(type[oldIndices], members[oldIndices])]
        end
        # calculate the new composite operator
        newOperator = prod([operators[(string(type[i]), [members[i]])] for i in newIndices])

        operators[(type, members)] = oldOperator * newOperator
    end
    return operators
end
export CreateProductOperator


"""
Given c^† matrix, construct c matrix, c^†c matrix and c c^† matrix.
"""
function CreateDNH(
        operators::Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}},
        site::Int64,
    )
    operators[("-", [site])] = operators[("+", [site])]'
    operators[("n", [site])] = operators[("+", [site])] * operators[("-", [site])]
    operators[("h", [site])] = operators[("-", [site])] * operators[("+", [site])]
    return operators
end
export CreateDNH


"""
Multiply individual requirements on occupancy and magnetisation
into a single if-condition.
"""
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
        requirement = nothing
    end
    return requirement
end
export CombineRequirements


"""
Given the symmetry requirements and a basis, return the set of quantum numbers
for each vector in the basis.
"""
function QuantumNosForBasis(
        currentSites::Vector{Int64}, 
        symmetries::Vector{Char},
        basisStates::Vector{Dict{BitVector, Float64}},
    )

    # obtain the symmetry operators as requested.
    symmetryOperators = Vector{Tuple{String,Vector{Int64},Float64}}[]

    # if 'N' is provided, total number operator commutes.
    if 'N' in symmetries
        push!(symmetryOperators, [("n", [i], 1.) for i in eachindex(currentSites)])
    end
    # if 'S' is provided, total Sz operator commutes.
    if 'S' in symmetries
        push!(symmetryOperators, [("n", [i], (-1)^(i+1)) for i in eachindex(currentSites)])
    end

    if !isempty(symmetries)
        
        # for each vector |Ψ⟩ in the basis, calculate the quantum numbers using
        # (⟨Ψ|O_1|Ψ⟩, ⟨Ψ|O_2|Ψ⟩), where O_1 and O_2 are the symmetry operators.
        return [Tuple(round(Int, GenCorrelation(state, operator)) for operator in symmetryOperators) 
                      for state in basisStates]
    else
        return nothing
    end
end
export QuantumNosForBasis


"""
Obtain operator matrices for the fundamental operators like c, c^†, n etc.
"""
function InitiateMatrices(currentSites, hamltFlow, initBasis)

    # matrix for the bond antisymmetriser operator, that is used to "attach" new sites
    bondAntiSymmzer = length(currentSites) == 1 ? sigmaz : kron(fill(sigmaz, length(currentSites))...)

    # zero hamiltonian with the appropriate size
    hamltMatrix = diagm(fill(0.0, 2^length(currentSites)))

    # create the c and cdaggers
    create, basket, newSitesFlow = UpdateRequirements(hamltFlow)
    operators = Dict{Tuple{String, Vector{Int64}}, Matrix{Float64}}()
    for site in currentSites 
        operators[("+", [site])] = OperatorMatrix(initBasis, [("+", [site], 1.0)])
        operators = CreateDNH(operators, site)
    end
    return operators, bondAntiSymmzer, hamltMatrix, newSitesFlow, create, basket
end
export InitiateMatrices


"""
Create filenames for saving data. Also make the initial write.
"""
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
export SetupDataWrite


"""
Obtain spectrum of Hamiltonian, making use of symmetries if available.
"""
function Diagonalise(
        hamltMatrix::Matrix{Float64},
        quantumNos::Union{Nothing, Vector{NTuple{1, Int64}}, Vector{NTuple{2, Int64}}},
    )
    # initialise matrix and vector to store eigenvectors and eigenvalues
    eigenVecs = zeros(size(hamltMatrix)...)
    eigVals = zeros(size(hamltMatrix)[2])

    if !isnothing(quantumNos)

        # loop over symmetry sectors
        for quantumNo in unique(quantumNos)

            # obtain the states for this sector
            indices = findall(==(quantumNo), quantumNos)

            # diagonalise this smaller subHamiltonian
            F = eigen(Hermitian(hamltMatrix[indices, indices]))
            eigenVecs[indices, indices] .= F.vectors
            eigVals[indices] .= F.values
        end

        # having consolidated the spectrum from all sectors, sort them
        # according to the energy eigenvalues.
        sortSequence = sortperm(eigVals)
        eigenVecs = eigenVecs[:, sortperm(eigVals)]
        quantumNos = quantumNos[sortperm(eigVals)]
        eigVals = sort(eigVals)
    else

        # if no symmetries, just diagonalise the whole matrix
        F = eigen(Hermitian(hamltMatrix))
        eigenVecs .= F.vectors
        eigVals .= F.values
    end
    return eigVals, eigenVecs, quantumNos
end
export Diagonalise


"""
Retain a given number of low energy states from the spectrum, in order to
keep the size manageable.
"""
function TruncateSpectrum(
        corrQuantumNoReq::Union{Nothing,Function},
        rotation::Matrix{Float64},
        eigVals::Vector{Float64},
        maxSize::Int64,
        degenTol::Float64,
        currentSites::Vector{Int64},
        quantumNos::Union{Nothing,Vector{NTuple{1, Int64}},Vector{NTuple{2, Int64}}}, 
        maxMaxSize::Int64,
    )
    if isnothing(corrQuantumNoReq)
        
        retainIndices = 1:length(eigVals)

        # ensure we stop not exactly at maxSize, but at the last state
        # with the same energy as that at maxSize. This ensures we are
        # not truncating in the middle of a degeneracy, at the cost of
        # a slower next step.
        filter!(i -> eigVals[i] ≤ eigVals[maxSize] + degenTol, retainIndices)
    else

        # if only certain symmetries are needed, keep lowest lying states
        # only from those sectors.
        retainIndices = findall(q -> corrQuantumNoReq(q, maximum(currentSites)), quantumNos)
        otherIndices = findall(q -> !corrQuantumNoReq(q, maximum(currentSites)), quantumNos)

        # if we can afford to keep more states even after taking all states
        # from  the specified sectors, add the lowest-lying balance number 
        # of states from the other sectors.
        if length(retainIndices) < maxSize
            maxSizeRemain = maxSize - length(retainIndices)
            cutOffEnergy = eigVals[otherIndices][maxSizeRemain] + degenTol
            filter!(i -> eigVals[i] ≤ cutOffEnergy, otherIndices)

            append!(retainIndices, otherIndices)
        else
            cutOffEnergy = eigVals[retainIndices][maxSize] + degenTol
            filter!(i -> eigVals[i] ≤ cutOffEnergy, retainIndices)
        end
    end

    # ensure that the increased number of states due to degeneracy
    # does not cross the absolute upper cutoff maxMaxSize
    if length(retainIndices) > maxMaxSize && maxMaxSize > 0
        retainIndices = retainIndices[1:maxMaxSize]
    end

    rotation = rotation[:, retainIndices]
    if !isnothing(quantumNos)
        quantumNos = quantumNos[retainIndices]
    end
    eigVals = eigVals[retainIndices]
    return rotation, eigVals, quantumNos
end
export TruncateSpectrum


"""
Rotate and enlarge old operators after each diagonalisation.
"""
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

    # remove operators that no longer need to be tracked
    for k in setdiff(keys(operators), newBasket)
        delete!(operators, k)
    end

    # rotate and enlarge qubit operators of current system
    for key in keys(operators)
        operators[key] = kron(rotation' * operators[key] * rotation, identityEnv)
    end
    bondAntiSymmzer = rotation' * bondAntiSymmzer * rotation

    # rotate and enlarge operator needed to compute correlations
    for (name, corrOperator) in collect(corrOperatorDict)
        if !isnothing(corrOperator)
            corrOperatorDict[name] = kron(rotation' * corrOperator * rotation, identityEnv)
        end
    end
    return hamltMatrix, operators, bondAntiSymmzer, corrOperatorDict
end
export UpdateOldOperators


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
    maxMaxSize::Int64,
)

    # ensure each term of the Hamiltonian is sorted in indices
    for (i, hamlt) in enumerate(hamltFlow)
        for (j, (operatorString, members, coupling)) in enumerate(hamlt)
            
            newString, newMembers, sign = OrganiseOperator(operatorString, members)
            hamltFlow[i][j] = (newString, newMembers, coupling * sign)
        end
    end

    # get the sites present at the zeroth step. use that to 
    # construct the initial basis.
    currentSites = collect(1:maximum(maximum.([opMembers for (_, opMembers, _) in hamltFlow[1]])))
    initBasis = BasisStates(maximum(currentSites))

    # get the quantum numbers for each state in the basis
    quantumNos = QuantumNosForBasis(currentSites, symmetries, initBasis)

    # get the filenames for saving data
    savePaths = SetupDataWrite(dataDir, hamltFlow, initBasis, currentSites, symmetries)

    # get the initial matrices, the sequence of new sites
    # being added at each step, and the various operators
    # that must be created at each step.
    operators, bondAntiSymmzer, hamltMatrix, newSitesFlow, create, basket = InitiateMatrices(currentSites, hamltFlow, initBasis)

    # if there are correlations, add them to the set of
    # operators that must be created and updated
    for (name, correlationDef) in correlationDefDict
        correlationMaxMember = maximum([maximum(members) for (_, members, _) in correlationDef])
        corrDefFlow = [ifelse(correlationMaxMember ∈ newSites, correlationDef, eltype(correlationDef)[]) for newSites in newSitesFlow]

        corrCreate, corrBasket = UpdateRequirements(corrDefFlow, newSitesFlow)
        create = [[c1; c2] for (c1, c2) in zip(create, corrCreate)]
        basket = [[c1; c2] for (c1, c2) in zip(basket, corrBasket)]
    end

    # create the first batch of operators
    for operator in create[1]
        operators = CreateProductOperator(operator, operators)
    end

    # stores operators that are necessary for calculating correlations
    corrOperatorDict = Dict{String, Union{Nothing, Matrix{Float64}}}(name => nothing for name in keys(correlationDefDict))

    # if any correlation operator can be calculated right away, do it
    for (name, correlationDef) in correlationDefDict
        corrDefKeys = [(type, members) for (type, members, _) in correlationDef]
        if isnothing(corrOperatorDict[name]) && all(∈(keys(operators)), corrDefKeys)
            corrOperatorDict[name] = sum([coupling * operators[(type, members)] for (type, members, coupling) in correlationDef])
        end
    end

    specFuncOperators = nothing
    if !isempty(specFuncNames)
        specFuncOperators = Dict{String, Vector{Matrix{Float64}}}(name => Matrix{Float64}[] for name in specFuncNames)
    end

    # if any spectral function operators can be calculated right away, do it
    if !isempty(specFuncNames) && all(name -> !isnothing(corrOperatorDict[name]), specFuncNames)
        for name in filter(name -> !isnothing(corrOperatorDict[name]), specFuncNames)
            push!(specFuncOperators[name], corrOperatorDict[name])
        end
    end

    resultsDict = Dict{String, Union{Nothing, Float64}}(name => nothing for name in keys(correlationDefDict))
    @assert "energyPerSite" ∉ keys(resultsDict)
    resultsDict["energyPerSite"] = nothing

    pbar = Progress(length(hamltFlow); enabled=!silent)
    for (step, hamlt) in enumerate(hamltFlow)
        # construct the Hamiltonian for this step
        for (type, members, strength) in hamlt
            hamltMatrix += strength * operators[(type, members)]
        end
        
        # get spectrum
        eigVals, rotation, quantumNos = Diagonalise(hamltMatrix, quantumNos)

        if step == length(hamltFlow)
            # if this is the last step, save data and calculate any correlations
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
            if !isnothing(corrQuantumNoReq) && !isnothing(quantumNos)
                indices = findall(q -> corrQuantumNoReq(q, maximum(currentSites)), quantumNos)
            end
            finalState = rotation[:, indices[sortperm(eigVals[indices])[1]]]

            # calculate correlations using the ground state of the final step
            for (name, correlationDef) in correlationDefDict
                if !isnothing(corrOperatorDict[name])
                    resultsDict[name] = finalState' * corrOperatorDict[name] * finalState
                end
            end

            next!(pbar; showvalues=[("Size", size(hamltMatrix))])
            break
        end

        # construct a basis and identity matrix for the new sites
        newBasis = BasisStates(length(newSitesFlow[step+1]))

        identityEnv = length(newSitesFlow[step+1]) == 1 ? I(2) : kron(fill(I(2), length(newSitesFlow[step+1]))...)

        # save data
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

        # truncate the spectrum
        if length(eigVals) > div(maxSize, 2^length(newSitesFlow[step+1]))
            rotation, eigVals, quantumNos = TruncateSpectrum(quantumNoReq, rotation, eigVals, 
                                                             div(maxSize, 2^length(newSitesFlow[step+1])), 
                                                             degenTol, currentSites, quantumNos, 
                                                             div(maxMaxSize, 2^length(newSitesFlow[step+1])),
                                                            )
        end

        # update the quantum numbers, due to newly added sites
        if !isnothing(quantumNos)
            newQuantumNos = QuantumNosForBasis(newSitesFlow[step+1], symmetries, newBasis)

            quantumNos = vcat([[quantumNo .+ newQuantumNo for newQuantumNo in newQuantumNos]
                           for quantumNo in quantumNos]...)
        end

        # rotate and enlarge existing operators
        hamltMatrix, operators, bondAntiSymmzer, corrOperator = UpdateOldOperators(eigVals, identityEnv, basket[step+1], 
                                                                                   operators, rotation, bondAntiSymmzer, corrOperatorDict
                                                                                  )

        # define the qbit operators for the new sites
        for site in newSitesFlow[step+1] 
            operators[("+", [site])] = kron(bondAntiSymmzer, OperatorMatrix(newBasis, [("+", [site - length(currentSites)], 1.0)]))
            operators = CreateDNH(operators, site)
        end

        # create new product operators that involve the new sites
        for operator in create[step+1]
            operators = CreateProductOperator(operator, operators; newSites=newSitesFlow[step+1])
        end

        # construct any new correlation operators that became available at this step
        for (name, correlationDef) in correlationDefDict
            corrDefKeys = [(type, members) for (type, members, _) in correlationDef]
            if isnothing(corrOperatorDict[name]) && all(∈(keys(operators)), corrDefKeys)
                corrOperatorDict[name] = sum([coupling * operators[(type, members)] for (type, members, coupling) in correlationDef])
            end
        end

        # if new spectral function operators can be created now, do it
        if !isempty(specFuncNames) && all(name -> !isnothing(corrOperatorDict[name]), specFuncNames)
            for name in specFuncNames
                push!(specFuncOperators[name], corrOperatorDict[name])
            end
        end

        # expand the antisymmetrizer
        bondAntiSymmzer = kron(bondAntiSymmzer, fill(sigmaz, length(newSitesFlow[step+1]))...)

        # update the current sites by adding the new sites
        append!(currentSites, newSitesFlow[step+1])

        # udpate progressbar
        next!(pbar; showvalues=[("Size", size(hamltMatrix))])
        GC.gc()
    end
    
    return savePaths, resultsDict, specFuncOperators
end
export IterDiag


function IterDiag(
    hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
    maxSize::Int64;
    occReq::Union{Nothing,Function}=nothing,
    # (o,N) -> div(N,2)-2 ≤ o ≤ div(N,2)+2 ## close to half-filling
    magzReq::Union{Nothing,Function}=nothing,
    # (m,N) -> m == N ## maximally polarised states
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
    maxMaxSize::Int64=0,
)
    @assert maxMaxSize == 0 || maxMaxSize ≥ maxSize

    exitCode = 0

    # get a list of the names of all the quantities that
    # have been requested. This involves ground state energy,
    # correlations, VNE, I2 and spec func.
    retainKeys = ["energyPerSite"]
    for quant in [correlationDefDict, vneDefDict, mutInfoDefDict, specFuncDefDict]
        append!(retainKeys, [k for k in keys(quant)])
    end

    # ensure that all the keys are unique,
    # otherwise one will be overwritten by another.
    @assert allunique(retainKeys)

    # ensure that there's at least two Hamiltonians to iterate over
    @assert length(hamltFlow) > 1

    # ensure that the symmetries are of the correct form (N and/or S),
    # and that the symmetry requirements are consistent with that
    # (the assert fails, for example, if we specify the occupancy
    # but don't pass 'N' as a symmetry)
    @assert all(∈("NS"), symmetries)
    if !isnothing(occReq)
        @assert 'N' in symmetries
    end
    if !isnothing(magzReq)
        @assert 'S' in symmetries
    end

    # reduce the mutual information calculation to that of VNE.
    # That is, if a certain I2(A:B) is requested, we convert that
    # to a request for S(A), S(B) and S(A:B). This prevents us
    # from having to deal separately with I2 during the iteration.
    # We just need to keep track of which three S goes into a specific I2
    mutInfoToVneMap = Dict{String, Union{Nothing,NTuple{3, String}}}(name => nothing for name in keys(mutInfoDefDict))
    for (name, (sysA, sysB)) in mutInfoDefDict

        # for each mut info request, udpate the vne requirements
        subsystemVneNames = []

        # loop over the three VNE that we must compute 
        # in order to compute this I2: S(A), S(B), S(A ∪ B)
        for subsystem in [sysA, sysB, vcat(sysA, sysB)]

            if subsystem ∈ values(vneDefDict)
                # if the requested VNE has also been separately
                # requested by the user, no need for extra computations,
                # just add the name to the records for this I2
                subsystemVneName = filter(k -> vneDefDict[k] == subsystem, collect(keys(vneDefDict)))[1]
            else 
                # otherwise, add this new VNE calculation to
                # the list of VNE that must be computed
                subsystemVneName = randstring(5)
                while subsystemVneName ∈ keys(vneDefDict)
                    subsystemVneName = randstring(5)
                end
                vneDefDict[subsystemVneName] = subsystem
            end
            push!(subsystemVneNames, subsystemVneName)
        end
        mutInfoToVneMap[name] = Tuple(subsystemVneNames)
    end

    # just like I2 calculations were reduced to VNE calculations,
    # we will now reduce VNE calculations to just correlation calculations;
    # this means we only need to track correlation functions during
    # the iterative diagonalisation, makes things simple. The way to
    # reduce a VNE calculation to a correlation function calculation
    # is by noting that the matrix elements of the reduced density matrix
    # can be written as ρ_A(i,j) = ⟨Ψ| |A_i⟩⟨A_j| |Ψ>, where |A_i⟩ are the
    # basis vectors for the subspace A. The projectors/transition operators
    # |A_i⟩⟨A_j| are simply various combinations of c^† c. For example,
    # |0⟩⟨1| = c, |1⟩⟨1|=c^†c, etc. Therefore, each matrix element is a
    # static correlation corresponding to these operators.
    vneOperatorToCorrMap = Dict{String, Dict{Tuple{String, Vector{Int64}}, Union{Nothing,String}}}(name => Dict() for name in keys(vneDefDict))
    for (name, sites) in vneDefDict

        # create all the transition operators |A_i⟩⟨A_j| associated with
        # this subspace. Add these correlation requirements to the list.
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

    # similarly, add the spectral function operators to the correlation
    # list, this enables that these operators will get created and
    # udpated during the iterative diagonalisation. We will finally
    # not use these operators to compute a ground state correlation,
    # but instead employ them for spectral function calculations.
    specFuncToCorrMap = Dict{String, String}()
    for (name, operator) in specFuncDefDict
        specFuncToCorrMap[name] = randstring(5)
        correlationDefDict[specFuncToCorrMap[name]] = operator
    end

    # once we have obtained the complete list of correlation operators
    # that must be computed (From various sources such as entanglement
    # and spectral function), we sort the indices of these operators.
    for (k, v) in deepcopy(correlationDefDict)
        for (j, (operatorString, members, coupling)) in enumerate(v)
            newString, newMembers, sign = OrganiseOperator(operatorString, members)#, newSitesFlow[i])
            correlationDefDict[k][j] = (newString, newMembers, coupling * sign)
        end
    end

    # consolidate the quantum number requirements (such as
    # a specific filling or a specific magnetisation sector)
    # into a single if-condition.
    quantumNoReq = CombineRequirements(occReq, magzReq)
    corrQuantumNoReq = CombineRequirements(corrOccReq, corrMagzReq)

    # maps transition operators (such as (1, 0) ≡ |1⟩⟨0|) 
    # to second quantised operators
    operatorMap = Dict(
                       (1,1) => "n",
                       (0,0) => "h",
                       (1,0) => "+",
                       (0,1) => "-",
                      )


    # perform the actual iterative diagonalisation. savePaths
    # contains the filepaths where data is saved; resultsDict
    # contains energyPerSite and correlation values;
    # specFuncOperators contains the updated forms of the
    # operators that will be used for spectral function calculations.
    savePaths, resultsDict, specFuncOperators = IterDiag(hamltFlow, maxSize, symmetries,
                                                         correlationDefDict,
                                                         quantumNoReq, 
                                                         corrQuantumNoReq,
                                                         degenTol,
                                                         dataDir,
                                                         silent,
                                                         collect(values(specFuncToCorrMap)),
                                                         maxMaxSize,
                                                        )

    # if specFuncOperators was requested, convert them back to 
    # their original names, because they had been passed into 
    # the iterative diagonaliser with random names in order to
    # avoid overwriting.
    if !isempty(specFuncToCorrMap)
        specFuncOperators = Dict{String,Vector{Matrix{Float64}}}(name => specFuncOperators[corrName] for (name, corrName) in specFuncToCorrMap)
    end

    # if vne was requested, we now have the various matrix elements
    # of the RDM; all that's left is to reconstruct the RDM. For this,
    # we loop over all basis vectors of the reduced subspace, obtain
    # the transition operator associated with matrix element, track
    # track the corresponding correlation function value, and set that
    # as the RDM element at that position.
    for (name, sites) in vneDefDict
        reducedDM = zeros(2^length(sites), 2^length(sites))
        basisVectors = collect(Iterators.product(repeat([[1, 0]], length(sites))...))
        for (i, sequence_i) in enumerate(basisVectors)
            for j in eachindex(basisVectors)[i:end]
                sequence_j = basisVectors[j]

                # sequence_i and sequence_j are the incoming and outgoing basis vectors
                # of the transition. For a two qubit subspace, this transition can, for
                # example, be |(1, 1)⟩⟨(0, 0)| = |1⟩⟨0| × |1⟩⟨0| . In second quantisation, 
                # this is c^† c^†, so the `operatorChars` variable is "++" ("+" is the symbol
                # for creation operator in our formulation).
                operatorChars = prod([operatorMap[(s_i, s_j)] for (s_i, s_j) in zip(sequence_i, sequence_j)])

                # this returns the specific name given to the computation request
                # for this particular operator string,during the iterative diag.
                corrName = vneOperatorToCorrMap[name][(operatorChars, sites)]

                # having obtained the correlation name, we can easily obtain its value
                # by getting the corresponding key from resultsDict.
                matrixElement = resultsDict[corrName]
                reducedDM[i, j] = matrixElement
                reducedDM[j, i] = matrixElement'
            end
        end

        # if trace of RDM is zero, something is wrong,
        # set exit code to non-zero to flag it,
        # otherwise normalise it.
        if tr(reducedDM) < 1e-10
            exitCode = 1
        else
            reducedDM /= tr(reducedDM)
        end

        # diagonalise RDM and use non-zero eigenvalues
        # to calculate VNE
        rdmEigVals = filter(>(0), eigvals(Hermitian(reducedDM)))
        resultsDict[name] = -sum(rdmEigVals .* log.(rdmEigVals))
    end

    # in the block above, we have already all VNE, including those
    # necessary for I2. In order to calculate I2, we simply track
    # the codenames of the VNE required for a specific I2, and 
    # combine them appropriately.
    for name in keys(mutInfoDefDict)
        vneNames = mutInfoToVneMap[name]
        resultsDict[name] = resultsDict[vneNames[1]] + resultsDict[vneNames[2]] - resultsDict[vneNames[3]]
        if resultsDict[name] < -1e-10
            exitCode = 2
        end
    end

    # delete all keys that were not specifically requested.
    # this deletes intermediate quantities that might have
    # been generated in the process.
    for name in keys(resultsDict)
        if name ∉ retainKeys
            delete!(resultsDict, name)
        end
    end

    # if there's entanglement involved, return exit code. if spectral function,
    # return that. Otherwise just return the results.
    if isempty(vneDefDict) && isempty(mutInfoDefDict) && isnothing(specFuncOperators)
        return savePaths, resultsDict
    elseif (!isempty(vneDefDict) || !isempty(mutInfoDefDict)) && isnothing(specFuncOperators)
        return savePaths, resultsDict, exitCode
    elseif isempty(vneDefDict) && isempty(mutInfoDefDict) && !isnothing(specFuncOperators)
        return savePaths, resultsDict, specFuncOperators
    elseif (!isempty(vneDefDict) || !isempty(mutInfoDefDict)) && !isnothing(specFuncOperators)
        return savePaths, resultsDict, exitCode, specFuncOperators
    end
end
export IterDiag


function IterSpecFunc(
        savePaths::Vector{String},
        specFuncOperators::Dict{String, Vector{Matrix{Float64}}},
        freqValues::Vector{Float64},
        standDev::Union{Vector{Float64}, Float64};
        degenTol::Float64=0.,
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
        eigVecs = [collect(col) for col in eachcol(data["basis"])]
        eigVals = data["eigVals"]
        quantumNos = data["quantumNos"]
        currentSites = data["currentSites"]

        if isnothing(quantumNos)
            minimalEigVecs = eigVecs
            minimalEigVals = eigVals
        else
            lowEnergyIndices = findall(q -> quantumNoReq(q, length(currentSites)), quantumNos)
            groundStateEnergy = minimum(eigVals[lowEnergyIndices])
            filter!(i -> eigVals[i] ≤ groundStateEnergy + degenTol, lowEnergyIndices)
            excitedIndices = findall(q -> excQuantumNoReq(q, length(currentSites)), quantumNos)
            allowedIndices = [lowEnergyIndices; excitedIndices]

            minimalEigVecs = eigVecs[allowedIndices]
            minimalEigVals = eigVals[allowedIndices]
            @assert abs(groundStateEnergy - minimum(minimalEigVals)) < 1e-10
        end

        specFunc = SpecFunc(minimalEigVals, minimalEigVecs, 
                            Dict(name => specFuncOperators[name][i] for name in keys(specFuncOperators)),
                            freqValues, standDev;
                            normalise=true, degenTol=degenTol)
        totalSpecFunc .+= specFunc
    end
    return totalSpecFunc
end


"""
A major component of the iterative diagonalisation algorithm is the
need to create and update operators carefully. For example, if a term
in the Hamiltonian is c^†_1 c_2 c^†_3 c_4, with sites 1, 2 entering 
at the first step and 3, 4 entering at second step, this operator must
be created as follows: c^†_1 c_2 must be created at the beginning of the
first step and updated(=rotated) at the end of the first step, and 
c^†_3 c_4 must be created at the end of the first step (when 3 and 4
enter the game). Once these two operators are available, the total operator
is created through a tensor product. This becomes more complicated if the
total operator stretches across multiple steps; it leads to a graph like
structure where various parts of the operator need to be created and
updated at various steps. This functions takes care of that bookkeeping;
it calculates exactly which operators need to be created at every step,
and which operators must be kept updating for future use. If any intermediate
stage operator has served its purpose, it is dropped for optimal usage of
resources.
"""
function UpdateRequirements(
        hamltFlow::Vector{Vector{Tuple{String,Vector{Int64},Float64}}},
        newSitesFlow::Vector{Vector{Int64}},
    )

    # toRetain tracks all operators that must be kept "fresh" at any given step
    toRetain = Vector{Tuple{String,Vector{Int64}}}[[] for _ in hamltFlow]

    # toCreate tracks all operators that must be created at any given step
    toCreate = Vector{Tuple{String,Vector{Int64}}}[[] for _ in hamltFlow]

    # start from the last step of the Hamiltonian flow
    for step in reverse(eachindex(hamltFlow))
        newSites = newSitesFlow[step]
        
        # the new operators that appear at any step must of course be created 
        # at that step, so they go into the toCreate list.
        if !isempty(hamltFlow[step])
            append!(toCreate[step], [(type, members) for (type, members, _) in hamltFlow[step]])
        end

        # if I am at any step apart from the last, I need to
        # consider updating operators as well. At the last step, there's no
        # need to update operators because there are no furthe steps.
        if step < length(hamltFlow)

            # if a certain operator is designated to be retained
            # at the next step, then it must be either retained
            # or created at the present step.
            for (type, members) in toRetain[step+1]

                # if this operator involves the new sites added
                # at this step, then it must actually be created
                # at this step, hence add to toCreate list.
                if !isempty(intersect(members, newSites))
                    push!(toCreate[step], (type, members))
                else
                    # otherwise just retain
                    push!(toRetain[step], (type, members))
                end
            end
        end

        # for the operators that must be created at the current
        # step, we must ensure that its "old components" are
        # created/retained at previous steps.
        for (type, members) in toCreate[step]
            oldMembers = findall(∉(newSites), members)

            # if the operator has old sites, then extract the
            # operator for the old sites and add to toRetain list,
            # since the old composite operator must be kept fresh
            # for multiplying with the new composite operator in
            # order to give the total operator.
            if !isempty(oldMembers)
                @assert oldMembers[end] == oldMembers[1] + length(oldMembers) - 1
                push!(toRetain[step], (type[oldMembers], members[oldMembers]))
            end
        end

        # drop any duplicates that might have been added
        toRetain[step] = unique(toRetain[step])
    end
    return toCreate, toRetain, newSitesFlow
end
export UpdateRequirements


"""
Extension of the previous function to the case when newSitesFlow
is embedded inside the Hamiltonian flow. This first reconstructs
the new sites being added at each step, then calls the previous
function.
"""
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


"""
Another extension of the previous function to the case when 
we want to determine which operators must be created and retaine
at every step in order to properly calculate a correlation function ⟨A⟩.
This requires constructing a dummy flow of the correlation operator A:
[0, 0, ... A, A, ...], where the first A appears when all sites occuring
in A have been added into the system.
"""
function UpdateRequirements(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        newSitesFlow::Vector{Vector{Int64}},
    )
    operatorMaxMember = maximum([maximum(members) for (_, members, _) in operator])
    operatorFlow = [ifelse(operatorMaxMember ∈ newSites, operator, Tuple{String,Vector{Int64},Float64}[]) for newSites in newSitesFlow]

    return UpdateRequirements(operatorFlow, newSitesFlow)
end
export UpdateRequirements


"""
Slice a Hamiltonian into multiple parts in order to
constitute a Hamiltonian flow that can be studied
using iterative diagonalisation. The slicing happens
at the indices provided in indexParitions; for each
index i, we retain all terms in the Hamiltonian that
act on i. If this is the first index, we retain all
terms that act only on indices upto i. For example, given
a Hamiltonian H = c^†_1 c_2 + c^†_2 c_3 + c^†_3 c_4 + h.c.
and indexPartitions = [2, 4], this function returns two
Hamiltonians, one for the first step H1 = c^†_1 c_2 and
one for the second step H2 = c^†_2 c_3 + c^†_3 c_4.
"""
function MinceHamiltonian(
        hamiltonian::Vector{Tuple{String, Vector{Int64}, Float64}},
        indexPartitions::Vector{Int64},
    )
    
    # first create the subspaces by defining start and stop
    # indices for each subspace.
    subspaces = NTuple{2, Int64}[]
    leadingIndex = 0
    for index in indexPartitions
        push!(subspaces, (leadingIndex+1, index))
        leadingIndex = index
    end

    # check that the subspaces are mutually exclusive
    @assert issorted(vcat(subspaces...))

    # allocate the Hamiltonian for each step of partitioning
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[[] for _ in indexPartitions]
    for (opType, members, coupling) in hamiltonian

        # for each term in the Hamiltonian, check in which subspace
        # its farthest index lies.
        subspaceIndex = findfirst(s -> s[1] ≤ maximum(members) ≤ s[2], subspaces)
        push!(hamFlow[subspaceIndex], (opType, members, coupling))
    end
    return hamFlow
end
export MinceHamiltonian
