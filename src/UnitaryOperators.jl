function generateOneLoopParticleTerms(IOMposition::Integer, alpha::Float64, num_entangled::Integer)
    operatorList = Dict{Tuple{String,Vector{Integer}},Float64}()
    for i in 1:num_entangled
        mergewith!(+, operatorList, Dict(("n+-", [1, 1 + 2 * i, IOMposition]) => 0.25 * alpha))
        mergewith!(+, operatorList, Dict(("n+-", [2, 1 + 2 * i, IOMposition]) => -0.25 * alpha))
        mergewith!(+, operatorList, Dict(("n+-", [1, 2 + 2 * i, IOMposition + 1]) => -0.25 * alpha))
        mergewith!(+, operatorList, Dict(("n+-", [2, 2 + 2 * i, IOMposition + 1]) => 0.25 * alpha))
        mergewith!(+, operatorList, Dict(("+-+-", [1, 2, 2 + 2 * i, IOMposition]) => 0.5 * alpha))
        mergewith!(+, operatorList, Dict(("+-+-", [2, 1, 1 + 2 * i, IOMposition + 1]) => 0.5 * alpha))
    end
    return operatorList
end


function generateOneLoopHoleTerms(IOMposition::Integer, alpha::Float64, num_entangled::Integer)
    operatorList = Dict{Tuple{String,Vector{Integer}},Float64}()
    for i in 1:num_entangled
        mergewith!(+, operatorList, Dict(("n+-", [1, IOMposition, 1 + 2 * i]) => 0.25 * alpha))
        mergewith!(+, operatorList, Dict(("n+-", [2, IOMposition, 1 + 2 * i]) => -0.25 * alpha))
        mergewith!(+, operatorList, Dict(("n+-", [1, IOMposition + 1, 2 + 2 * i]) => -0.25 * alpha))
        mergewith!(+, operatorList, Dict(("n+-", [2, IOMposition + 1, 2 + 2 * i]) => 0.25 * alpha))
        mergewith!(+, operatorList, Dict(("+-+-", [2, 1, IOMposition, 2 + 2 * i]) => 0.5 * alpha))
        mergewith!(+, operatorList, Dict(("+-+-", [1, 2, IOMposition + 1, 1 + 2 * i]) => 0.5 * alpha))
    end
    return operatorList
end


function KondoUnitaries(alpha::Float64, num_entangled::Integer, sectors::String)
    @assert sectors ∈ ("p", "h", "ph")
    IOMposition = 2 * num_entangled + 1
    if sectors == "p"
        return generateOneLoopParticleTerms(IOMposition, alpha, num_entangled)
    elseif sectors == "h"
        return generateOneLoopHoleTerms(IOMposition, alpha, num_entangled)
    else
        particleTerms = generateOneLoopParticleTerms(IOMposition, alpha, num_entangled)
        holeTerms = generateOneLoopParticleTerms(IOMposition + 2, alpha, num_entangled)
        return mergewith(+, particleTerms, holeTerms)
    end
    return operatorList
end
