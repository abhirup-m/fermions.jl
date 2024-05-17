function kondoKSpace(dispersionDict::Dict{Int64,Float64}, kondoDict::Dict{Tuple{Int64,Int64},Float64}, bathIntDict::Dict{Tuple{Int64,Int64,Int64,Int64},Float64})
    operatorList = Dict{Tuple{String,Vector{Int64}},Float64}()
    # for index in 0:3
    #     merge!(+, operatorList, Dict(("n", [2 * index + 1]) => TOLERANCE * 1000))
    #     merge!(+, operatorList, Dict(("n", [2 * index + 2]) => -(TOLERANCE * 1000)))
    # end

    for (momIndex, energy) in dispersionDict
        # if energy == 0
        #     continue
        # end
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
        # if abs(kondoIntVal) < TOLERANCE
        #     continue
        # end

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
        # if abs(bathIntVal) < TOLERANCE
        #     continue
        # end

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
    return operatorList
end


function kondoKSpace(dispersionDictArray::Vector{Dict{Int64,Float64}}, kondoDictArray::Vector{Dict{Tuple{Int64,Int64},Float64}}, bathIntDictArray::Vector{Dict{Tuple{Int64,Int64,Int64,Int64},Float64}})
    operatorListSet = Tuple{String,Vector{Int64}}[]
    couplingMatrix = Vector{Float64}[]
    for (dispersionDict, kondoDict, bathIntDict) in zip(dispersionDictArray, kondoDictArray, bathIntDictArray)
        operatorList = kondoKSpace(dispersionDict, kondoDict, bathIntDict)
        if isempty(operatorListSet)
            operatorListSet = collect(keys(operatorList))
        end
        push!(couplingMatrix, collect(values(operatorList)))
    end
    return operatorListSet, couplingMatrix
end
