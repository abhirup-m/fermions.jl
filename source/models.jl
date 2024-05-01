include("fermionise.jl")

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
    operatorList = Tuple{String,Float64,Vector{Int64}}[]
    mapSeq = Dict(index => i for (i, index) in enumerate(chosenIndices))

    dispersionDict = Dict(mapSeq[index] => dispersionArray[index] for index in chosenIndices)
    kondoCouplingDict = Dict(Tuple(mapSeq[p] for p in points) => kondoJArray[points...]
                             for points in Iterators.product(chosenIndices, chosenIndices))
    bathIntDict = Dict(Tuple(mapSeq[p] for p in points) => bathIntFunc(points)
                       for points in Iterators.product(chosenIndices, chosenIndices, chosenIndices, chosenIndices))
    # kinetic energy terms
    for momIndex in points
        push!(operatorList, ("n", dispersionDict[momIndex], [2 * momIndex + 1]))
        push!(operatorList, ("n", dispersionDict[momIndex], [2 * momIndex + 2]))
    end

    # Kondo terms, for all pairs of momenta
    for (momIndex1, momIndex2) in Iterators.product(points, points)

        # get the up and down indices for both momenta
        momUpIndex1, momUpIndex2 = 2 .* (momIndex1, momIndex2) .+ 1
        momDownIndex1, momDownIndex2 = 2 .* (momIndex1, momIndex2) .+ 2
        impUpIndex = 1
        impDownIndex = 2

        if kondoCouplingDict[(momIndex1, momIndex2)] == 0
            continue
        end

        # 1/4 J n_d↑ (c^†_{k1 ↑}c_{k2 ↑}) 
        push!(operatorList, ("n+-", 0.25 * kondoCouplingDict[(momIndex1, momIndex2)], [impUpIndex, momUpIndex1, momUpIndex2]))

        # 1/4 J n_d↑ (-c^†_{k1 ↓}c_{k2 ↓}) 
        push!(operatorList, ("n+-", -0.25 * kondoCouplingDict[(momIndex1, momIndex2)], [impUpIndex, momDownIndex1, momDownIndex2]))

        # -1/4 J n_d↓ (c^†_{k1 ↑}c_{k2 ↑}) 
        push!(operatorList, ("n+-", -0.25 * kondoCouplingDict[(momIndex1, momIndex2)], [impDownIndex, momUpIndex1, momUpIndex2]))

        # -1/4 J n_d↓ (-c^†_{k1 ↓}c_{k2 ↓}) 
        push!(operatorList, ("n+-", 0.25 * kondoCouplingDict[(momIndex1, momIndex2)], [impDownIndex, momDownIndex1, momDownIndex2]))

        # 1/2 c^†_{d ↑} c_{d ↓} c^†_{k1 ↓} c_{k2 ↑}
        push!(operatorList, ("+-+-", 0.5 * kondoCouplingDict[(momIndex1, momIndex2)], [impUpIndex, impDownIndex, momDownIndex1, momUpIndex2]))

        # 1/2 c^†_{d ↓} c_{d ↑} c^†_{k1 ↑} c_{k2 ↓}
        push!(operatorList, ("+-+-", 0.5 * kondoCouplingDict[(momIndex1, momIndex2)], [impDownIndex, impUpIndex, momUpIndex1, momDownIndex2]))
    end

    # bath interaction terms, for all quartets of momenta
    for momIndices in Iterators.product(points, points, points, points)

        # get the up and down indices for all momenta
        upIndices = collect(2 .* momIndices .+ 1)
        downIndices = collect(2 .* momIndices .+ 2)

        if bathIntDict[momIndices] == 0
            continue
        end

        # if length(Set(momIndices)) < 3
        #     continue
        # end

        # -0.5 W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↑} c_{4, ↑}
        push!(operatorList, ("+-+-", -0.5 * bathIntDict[momIndices], upIndices))

        # -0.5 W_{1,2,3,4} c^†_{1 ↓} c_{2, ↓} c^†_{3 ↓} c_{4, ↓}
        push!(operatorList, ("+-+-", -0.5 * bathIntDict[momIndices], downIndices))

        # W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↓} c_{4, ↓}
        push!(operatorList, ("+-+-", bathIntDict[momIndices], [upIndices[1:2]; downIndices[3:4]]))

    end

    return operatorList
end
