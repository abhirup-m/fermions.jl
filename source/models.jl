include("fermionise.jl")

function KondoKSpace(basisStates, numBathStates, couplings)
    # imup    imdn    k1up    k1dn    k2up    k2dn    k3up    k3dn ...
    # 0       1       2       3       4       5       6       7    ...

    if length(basisStates) == 0
        println("Chosen basis is empty, bailing out.")
        return
    end
    (dispersion_arr, kondoCouplingMatrix, bathIntFunc) = couplings
    @assert length(dispersion_arr) == numBathStates
    @assert size(kondoCouplingMatrix) == (numBathStates, numBathStates)

    operatorList = Tuple{String, Float64, Vector{Int64}}[]

    # kinetic energy terms
    for momIndex in 1:numBathStates
        push!(operatorList, ("n", dispersion_arr[momIndex], (2 * momIndex + 1)))
        push!(operatorList, ("n", dispersion_arr[momIndex], (2 * momIndex + 2)))
    end

    # Kondo terms, for all pairs of momenta
    for (momIndex1, momIndex2) in Iterators.product(1:numBathStates, 1:numBathStates)

        # get the up and down indices for both momenta
        momUpIndex1, momUpIndex2 = 2 .* (momIndex1, momIndex2) .+ 1
        momDownIndex1, momDownIndex2 = 2 .* (momIndex1, momIndex2) .+ 2
        impUpIndex = 1
        impDownIndex = 1

        # 1/4 J n_d↑ (c^†_{k1 ↑}c_{k2 ↑}) 
        push!(operatorList, ("n+-", 0.25 * kondoCouplingMatrix[momIndex1, momIndex2], [impUpIndex, momUpIndex1, momUpIndex2]))

        # 1/4 J n_d↑ (-c^†_{k1 ↓}c_{k2 ↓}) 
        push!(operatorList, ("n+-", -0.25 * kondoCouplingMatrix[momIndex1, momIndex2], [impUpIndex, momDownIndex1, momDownIndex2]))

        # -1/4 J n_d↓ (c^†_{k1 ↑}c_{k2 ↑}) 
        push!(operatorList, ("n+-", -0.25 * kondoCouplingMatrix[momIndex1, momIndex2], [impDownIndex, momUpIndex1, momUpIndex2]))

        # -1/4 J n_d↓ (-c^†_{k1 ↓}c_{k2 ↓}) 
        push!(operatorList, ("n+-", 0.25 * kondoCouplingMatrix[momIndex1, momIndex2], [impDownIndex, momDownIndex1, momDownIndex2]))

        # 1/2 c^†_{d ↑} c_{d ↓} c^†_{k1 ↓} c_{k2 ↑}
        push!(operatorList, ("+-+-", 0.5 * kondoCouplingMatrix[momIndex1, momIndex2], [impUpIndex, impDownIndex, momDownIndex1, momUpIndex2]))

        # 1/2 c^†_{d ↓} c_{d ↑} c^†_{k1 ↑} c_{k2 ↓}
        push!(operatorList, ("+-+-", 0.5 * kondoCouplingMatrix[momIndex1, momIndex2], [impDownIndex, impUpIndex, momUpIndex1, momDownIndex2]))
    end

    # bath interaction terms, for all quartets of momenta
    for momIndices in Iterators.product(1:numBathStates, 1:numBathStates, 1:numBathStates, 1:numBathStates)

        # get the up and down indices for all momenta
        upIndices = 2 .* momIndices .+ 1
        downIndices = 2 .* momIndices .+ 2

        # -0.5 W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↑} c_{4, ↑}
        push!(operatorList, ("+-+-", -0.5 * bathIntFunc[sort(collect(momIndices))], upIndices))

        # -0.5 W_{1,2,3,4} c^†_{1 ↓} c_{2, ↓} c^†_{3 ↓} c_{4, ↓}
        push!(operatorList, ("+-+-", -0.5 * bathIntFunc[sort(collect(momIndices))], downIndices))

        # W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↓} c_{4, ↓}
        push!(operatorList, ("+-+-", bathIntFunc[sort(collect(momIndices))], [upIndices[1:2]; downIndices[3:4]]))

    return generalOperatorMatrix(basisStates, operatorList)
end
