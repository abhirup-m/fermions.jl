function KondoUnitaries(alpha::Float64, num_entangled::Integer, sectors::String)
    @assert sectors ∈ ("p", "h", "ph")
    operatorList = Dict{Tuple{String,Vector{Int64}},Float64}()
    IOMposition = 2 * num_entangled + 1
    for sector in sectors
        if sector == 'h'
            for i in 1:num_entangled
                merge!(+, operatorList, Dict(("n+-", [1, IOMposition, 1 + 2 * i]) => 0.25 * alpha))
                merge!(+, operatorList, Dict(("n+-", [2, IOMposition, 1 + 2 * i]) => -0.25 * alpha))
                merge!(+, operatorList, Dict(("n+-", [1, IOMposition + 1, 2 + 2 * i]) => -0.25 * alpha))
                merge!(+, operatorList, Dict(("n+-", [2, IOMposition + 1, 2 + 2 * i]) => 0.25 * alpha))
                merge!(+, operatorList, Dict(("+-+-", [2, 1, IOMposition, 2 + 2 * i]) => 0.5 * alpha))
                merge!(+, operatorList, Dict(("+-+-", [1, 2, IOMposition + 1, 1 + 2 * i]) => 0.5 * alpha))
            end
        else
            for i in 1:num_entangled
                merge!(+, operatorList, Dict(("n+-", [1, 1 + 2 * i, IOMposition]) => 0.25 * alpha))
                merge!(+, operatorList, Dict(("n+-", [2, 1 + 2 * i, IOMposition]) => -0.25 * alpha))
                merge!(+, operatorList, Dict(("n+-", [1, 2 + 2 * i, IOMposition + 1]) => -0.25 * alpha))
                merge!(+, operatorList, Dict(("n+-", [2, 2 + 2 * i, IOMposition + 1]) => 0.25 * alpha))
                merge!(+, operatorList, Dict(("+-+-", [1, 2, 2 + 2 * i, IOMposition]) => 0.5 * alpha))
                merge!(+, operatorList, Dict(("+-+-", [2, 1, 1 + 2 * i, IOMposition + 1]) => 0.5 * alpha))
            end
        end
        IOMposition += 2
    end
        
    return operatorList
end
