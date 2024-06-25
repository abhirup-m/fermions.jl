function kondoKSpace(dispersionDict::Dict{Int64,Float64}, kondoDict::Dict{Tuple{Int64,Int64},Float64}, bathIntDict::Dict{Tuple{Int64,Int64,Int64,Int64},Float64}; bathField::Float64=0.0, tolerance::Float64=1e-16)
    operatorList = Dict{Tuple{String,Vector{Int64}},Float64}()

    for (momIndex, energy) in dispersionDict
        energy = abs(energy) > tolerance ? energy : 0
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 1]) => energy))
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 2]) => energy))
    end

    momenta_indices = 1:length(dispersionDict)
    indicesRepeatedTwice = ntuple(x -> momenta_indices, 2)
    indicesRepeatedFour = ntuple(x -> momenta_indices, 4)

    # Kondo terms, for all pairs of momenta
    for k_index_pair in Iterators.product(indicesRepeatedTwice...)

        # get the up and down indices for both momenta
        k_UpIndex1, k_UpIndex2 = 2 .* k_index_pair .+ 1
        k_DownIndex1, k_DownIndex2 = 2 .* k_index_pair .+ 2
        impUpIndex = 1
        impDownIndex = 2

        kondoIntVal = kondoDict[k_index_pair...]
        kondoIntVal = abs(kondoIntVal) > tolerance ? kondoIntVal : 0
        if kondoIntVal == 0
            continue
        end

        # 1/4 J n_d↑ (c^†_{k1 ↑}c_{k2 ↑}) 
        merge!(+, operatorList, Dict(("n+-", [impUpIndex, k_UpIndex1, k_UpIndex2]) => 0.25 * kondoIntVal))

        # 1/4 J n_d↑ (-c^†_{k1 ↓}c_{k2 ↓}) 
        merge!(+, operatorList, Dict(("n+-", [impUpIndex, k_DownIndex1, k_DownIndex2]) => -0.25 * kondoIntVal))

        # -1/4 J n_d↓ (c^†_{k1 ↑}c_{k2 ↑}) 
        merge!(+, operatorList, Dict(("n+-", [impDownIndex, k_UpIndex1, k_UpIndex2]) => -0.25 * kondoIntVal))

        # -1/4 J n_d↓ (-c^†_{k1 ↓}c_{k2 ↓}) 
        merge!(+, operatorList, Dict(("n+-", [impDownIndex, k_DownIndex1, k_DownIndex2]) => 0.25 * kondoIntVal))

        # 1/2 c^†_{d ↑} c_{d ↓} c^†_{k1 ↓} c_{k2 ↑}
        merge!(+, operatorList, Dict(("+-+-", [impUpIndex, impDownIndex, k_DownIndex1, k_UpIndex2]) => 0.5 * kondoIntVal))

        # 1/2 c^†_{d ↓} c_{d ↑} c^†_{k1 ↑} c_{k2 ↓}
        merge!(+, operatorList, Dict(("+-+-", [impDownIndex, impUpIndex, k_UpIndex1, k_DownIndex2]) => 0.5 * kondoIntVal))

    end

    # bath interaction terms, for all quartets of momenta
    for index4tuples in Iterators.product(indicesRepeatedFour...)
        bathIntVal = bathIntDict[index4tuples...]
        bathIntVal = abs(bathIntVal) > tolerance ? bathIntVal : 0
        if bathIntVal == 0
            continue
        end

        # get the up and down indices for all momenta
        upIndices = collect(2 .* index4tuples .+ 1)
        downIndices = collect(2 .* index4tuples .+ 2)

        # -0.5 W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↑} c_{4, ↑}
        merge!(+, operatorList, Dict(("+-+-", upIndices) => -0.5 * bathIntVal))

        # -0.5 W_{1,2,3,4} c^†_{1 ↓} c_{2, ↓} c^†_{3 ↓} c_{4, ↓}
        merge!(+, operatorList, Dict(("+-+-", downIndices) => -0.5 * bathIntVal))

        # W_{1,2,3,4} c^†_{1 ↑} c_{2, ↑} c^†_{3 ↓} c_{4, ↓}
        merge!(+, operatorList, Dict(("+-+-", [upIndices[1:2]; downIndices[3:4]]) => bathIntVal))

    end

    # conduction bath magnetic field
    for i in momenta_indices
        merge!(+, operatorList, Dict(("n", [2 * i + 1]) => bathField))
        merge!(+, operatorList, Dict(("n", [2 * i + 2]) => -bathField))
    end

    return operatorList
end


function kondoKSpace(dispersionDictArray::Vector{Dict{Int64,Float64}}, kondoDictArray::Vector{Dict{Tuple{Int64,Int64},Float64}}, bathIntDictArray::Vector{Dict{Tuple{Int64,Int64,Int64,Int64},Float64}}; bathField::Float64=0.0, tolerance::Float64=1e-16)
    operatorListSet = Tuple{String,Vector{Int64}}[]
    couplingMatrix = Vector{Float64}[]
    operatorListResults = fetch.([Threads.@spawn kondoKSpace(args...; bathField=bathField, tolerance=tolerance) for args in zip(dispersionDictArray, kondoDictArray, bathIntDictArray)])
    operatorListSet = collect(keys(operatorListResults[1]))
    couplingMatrix = [collect(values(operatorList)) for operatorList in operatorListResults]
    return operatorListSet, couplingMatrix
end


function kondoRSpace(dispersion::Vector{Float64}, kondoJ::Float64; tolerance::Float64=1e-16)
    operatorList = Dict{Tuple{String,Vector{Int64}},Float64}()
    for (momIndex, energy) in enumerate(dispersion)
        energy = abs(energy) > tolerance ? energy : 0
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 1]) => energy))
        merge!(+, operatorList, Dict(("n", [2 * momIndex + 2]) => energy))
    end
    momenta_indices = 1:length(dispersion)
    indicesRepeatedTwice = ntuple(x -> momenta_indices, 2)

    # Kondo terms, for all pairs of momenta
    for k_index_pair in Iterators.product(indicesRepeatedTwice...)

        # get the up and down indices for both momenta
        k_UpIndex1, k_UpIndex2 = 2 .* k_index_pair .+ 1
        k_DownIndex1, k_DownIndex2 = 2 .* k_index_pair .+ 2
        impUpIndex = 1
        impDownIndex = 2

        # 1/4 J n_d↑ (c^†_{k1 ↑}c_{k2 ↑}) 
        merge!(+, operatorList, Dict(("n+-", [impUpIndex, k_UpIndex1, k_UpIndex2]) => 0.25 * kondoJ))

        # 1/4 J n_d↑ (-c^†_{k1 ↓}c_{k2 ↓}) 
        merge!(+, operatorList, Dict(("n+-", [impUpIndex, k_DownIndex1, k_DownIndex2]) => -0.25 * kondoJ))

        # -1/4 J n_d↓ (c^†_{k1 ↑}c_{k2 ↑}) 
        merge!(+, operatorList, Dict(("n+-", [impDownIndex, k_UpIndex1, k_UpIndex2]) => -0.25 * kondoJ))

        # -1/4 J n_d↓ (-c^†_{k1 ↓}c_{k2 ↓}) 
        merge!(+, operatorList, Dict(("n+-", [impDownIndex, k_DownIndex1, k_DownIndex2]) => 0.25 * kondoJ))

        # 1/2 c^†_{d ↑} c_{d ↓} c^†_{k1 ↓} c_{k2 ↑}
        merge!(+, operatorList, Dict(("+-+-", [impUpIndex, impDownIndex, k_DownIndex1, k_UpIndex2]) => 0.5 * kondoJ))

        # 1/2 c^†_{d ↓} c_{d ↑} c^†_{k1 ↑} c_{k2 ↓}
        merge!(+, operatorList, Dict(("+-+-", [impDownIndex, impUpIndex, k_UpIndex1, k_DownIndex2]) => 0.5 * kondoJ))

    end
    return operatorList
end
