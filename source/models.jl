using ProgressMeter

function KondoKSpace(chosenIndices::Vector{Int64}, dispersionArray::Vector{Float64}, kondoJArray::Matrix{Float64}, bathIntFunc)
    # Hamiltonian is of the form
    # H = \sum_k ε_k n̂_k + \sum_{k1, k2} J_{k1, k2} S_d ⋅σ_{α,β} c^†_{k1 α} c_{k2 β}
    #     - 0.5\sum_{1,2,3,4} W_{1,2,3,4} \sum_σ (c^†_{1 σ} c_{2 σ}c^†_{3 σ} c_{4 σ} 
    #     - c^†_{1 σ} c_{2 σ}c^†_{3 σ̄} c_{4 σ̄})
    # Various sites and states are labelled as
    # imup    imdn    k1up    k1dn    k2up    k2dn    k3up    k3dn ...
    # 0       1       2       3       4       5       6       7    ...

    # We need the following couplings:
    # dispersion_arr: Dictionary of kinetic energies. They keys are integers i that index
    # the momentum states, and the values are the kinetic energies Ek[k_i] for those indices.
    # kondoCouplingMatrix: Dictionary that defines the Kondo interaction terms. The keys
    # are tuples of momentum indices (i, j), and the value at the key is the 
    # Kondo coupling J_{ki, kj}.
    # bathIntFunc: Dictionary that defines the 4-point bath interaction terms. The keys
    # are quartets of momentum indices [i,j,p,q], and the value at the key is the 
    # interaction strength W_{ki, kj, kp, kq}.

    @assert size(kondoJArray) == (length(dispersionArray), length(dispersionArray))

    points = 1:length(chosenIndices)
    operatorList = Dict{Tuple{String,Vector{Int64}},Float64}()
    mapSeq = Dict(index => i for (i, index) in enumerate(chosenIndices))

    dispersionDict = Dict(mapSeq[index] => dispersionArray[index] for index in chosenIndices)
    kondoCouplingDict = Dict(Tuple(mapSeq[p] for p in points) => kondoJArray[points...]
                             for points in Iterators.product(chosenIndices, chosenIndices))
    bathIntDict = Dict(Tuple(mapSeq[p] for p in points) => bathIntFunc(points)
                       for points in Iterators.product(chosenIndices, chosenIndices, chosenIndices, chosenIndices))
    # kinetic energy terms
    for momIndex in points
        kineticEnergy = round(dispersionDict[momIndex], digits=trunc(Int, -log10(TOLERANCE)))
        if kineticEnergy == 0
            continue
        end
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 1]) => kineticEnergy))
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 2]) => kineticEnergy))
    end

    # Kondo terms, for all pairs of momenta
    for (momIndex1, momIndex2) in Iterators.product(reverse(points), reverse(points))
        # get the up and down indices for both momenta
        momUpIndex1, momUpIndex2 = 2 .* (momIndex1, momIndex2) .+ 1
        momDownIndex1, momDownIndex2 = 2 .* (momIndex1, momIndex2) .+ 2
        impUpIndex = 1
        impDownIndex = 2

        kondoCouplingk1k2 = round(kondoCouplingDict[(momIndex1, momIndex2)], digits=trunc(Int, -log10(TOLERANCE)))
        if kondoCouplingk1k2 == 0
            continue
        end
        kondoCouplingk1k2 = 1

        # 1/4 J n_d↑ (c^†_{k1 ↑}c_{k2 ↑}) 
        # merge!(+, operatorList, Dict(("n+-", [impUpIndex, momUpIndex1, momUpIndex2]) => 0.25 * kondoCouplingk1k2))

        # 1/4 J n_d↑ (-c^†_{k1 ↓}c_{k2 ↓}) 
        # merge!(+, operatorList, Dict(("n+-", [impUpIndex, momDownIndex1, momDownIndex2]) => -0.25 * kondoCouplingk1k2))

        # -1/4 J n_d↓ (c^†_{k1 ↑}c_{k2 ↑}) 
        # merge!(+, operatorList, Dict(("n+-", [impDownIndex, momUpIndex1, momUpIndex2]) => -0.25 * kondoCouplingk1k2))

        # -1/4 J n_d↓ (-c^†_{k1 ↓}c_{k2 ↓}) 
        # merge!(+, operatorList, Dict(("n+-", [impDownIndex, momDownIndex1, momDownIndex2]) => 0.25 * kondoCouplingk1k2))

        # 1/2 c^†_{d ↑} c_{d ↓} c^†_{k1 ↓} c_{k2 ↑}
        merge!(+, operatorList, Dict(("+-+-", [impUpIndex, impDownIndex, momDownIndex1, momUpIndex2]) => 0.5 * kondoCouplingk1k2))

        # 1/2 c^†_{d ↓} c_{d ↑} c^†_{k1 ↑} c_{k2 ↓}
        merge!(+, operatorList, Dict(("+-+-", [impDownIndex, impUpIndex, momUpIndex1, momDownIndex2]) => 0.5 * kondoCouplingk1k2))

    end

    # bath interaction terms, for all quartets of momenta
    for momIndices in Iterators.product(points, points, points, points)
        bathIntCoup = round(bathIntDict[momIndices], digits=trunc(Int, -log10(TOLERANCE)))
        if bathIntCoup == 0
            continue
        end

        # get the up and down indices for all momenta
        upIndices = collect(2 .* momIndices .+ 1)
        downIndices = collect(2 .* momIndices .+ 2)

        # -0.5 W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↑} c_{4, ↑}
        merge!(+, operatorList, Dict(("+-+-", upIndices) => -0.5 * bathIntCoup))

        # -0.5 W_{1,2,3,4} c^†_{1 ↓} c_{2, ↓} c^†_{3 ↓} c_{4, ↓}
        merge!(+, operatorList, Dict(("+-+-", downIndices) => -0.5 * bathIntCoup))

        # W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↓} c_{4, ↓}
        merge!(+, operatorList, Dict(("+-+-", [upIndices[1:2]; downIndices[3:4]]) => bathIntCoup))

    end

    return operatorList
end


function KondoKSpace(chosenIndicesSet::Vector{Vector{Int64}}, dispersionArray::Vector{Float64}, kondoJArray::Matrix{Float64}, bathIntFunc)
    operatorListSet = Dict{Tuple{String,Vector{Int64}},Vector{Float64}}()
    for chosenIndices in chosenIndicesSet
        operatorList = KondoKSpace(chosenIndices, dispersionArray, kondoJArray, bathIntFunc)
        operatorListVectored = Dict(k => [v] for (k, v) in operatorList)
        mergewith!((v1, v2) -> [v1; v2], operatorListSet, operatorListVectored)
    end
    return operatorListSet
end


function kondoKSpace_test(chosenIndices, dispersionArray, kondoDict, bathIntDict)
    k_index_pairs = enumerate(chosenIndices)
    operatorList = Dict{Tuple{String,Vector{Int64}},Float64}()
    for (momIndex, momentum) in k_index_pairs
        if dispersionArray[momentum] == 0
            continue
        end
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 1]) => dispersionArray[momentum]))
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 2]) => dispersionArray[momentum]))
    end

    # Kondo terms, for all pairs of momenta
    for ((momIndex1, mom1), (momIndex2, mom2)) in Iterators.product(k_index_pairs, k_index_pairs)

        # get the up and down indices for both momenta
        momUpIndex1, momUpIndex2 = 2 .* (momIndex1, momIndex2) .+ 1
        momDownIndex1, momDownIndex2 = 2 .* (momIndex1, momIndex2) .+ 2
        impUpIndex = 1
        impDownIndex = 2

        kondok1k2 = kondoDict[(mom1, mom2)]
        if kondok1k2 == 0
            continue
        end

        # 1/4 J n_d↑ (c^†_{k1 ↑}c_{k2 ↑}) 
        merge!(+, operatorList, Dict(("n+-", [impUpIndex, momUpIndex1, momUpIndex2]) => 0.25 * kondok1k2))

        # 1/4 J n_d↑ (-c^†_{k1 ↓}c_{k2 ↓}) 
        merge!(+, operatorList, Dict(("n+-", [impUpIndex, momDownIndex1, momDownIndex2]) => -0.25 * kondok1k2))

        # -1/4 J n_d↓ (c^†_{k1 ↑}c_{k2 ↑}) 
        merge!(+, operatorList, Dict(("n+-", [impDownIndex, momUpIndex1, momUpIndex2]) => -0.25 * kondok1k2))

        # -1/4 J n_d↓ (-c^†_{k1 ↓}c_{k2 ↓}) 
        merge!(+, operatorList, Dict(("n+-", [impDownIndex, momDownIndex1, momDownIndex2]) => 0.25 * kondok1k2))

        # 1/2 c^†_{d ↑} c_{d ↓} c^†_{k1 ↓} c_{k2 ↑}
        merge!(+, operatorList, Dict(("+-+-", [impUpIndex, impDownIndex, momDownIndex1, momUpIndex2]) => 0.5 * kondok1k2))

        # 1/2 c^†_{d ↓} c_{d ↑} c^†_{k1 ↑} c_{k2 ↓}
        merge!(+, operatorList, Dict(("+-+-", [impDownIndex, impUpIndex, momUpIndex1, momDownIndex2]) => 0.5 * kondok1k2))

    end

    # bath interaction terms, for all quartets of momenta
    for pairs in Iterators.product(k_index_pairs, k_index_pairs, k_index_pairs, k_index_pairs)
        indices = [p[1] for p in pairs]
        momenta = tuple([p[2] for p in pairs]...)
        if abs(bathIntDict[momenta]) < TOLERANCE
            continue
        end

        # get the up and down indices for all momenta
        upIndices = collect(2 .* indices .+ 1)
        downIndices = collect(2 .* indices .+ 2)

        # -0.5 W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↑} c_{4, ↑}
        merge!(+, operatorList, Dict(("+-+-", upIndices) => -0.5 * bathIntDict[momenta]))

        # -0.5 W_{1,2,3,4} c^†_{1 ↓} c_{2, ↓} c^†_{3 ↓} c_{4, ↓}
        merge!(+, operatorList, Dict(("+-+-", downIndices) => -0.5 * bathIntDict[momenta]))

        # W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↓} c_{4, ↓}
        merge!(+, operatorList, Dict(("+-+-", [upIndices[1:2]; downIndices[3:4]]) => bathIntDict[momenta]))

    end
    return operatorList
end
