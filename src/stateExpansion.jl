function stateExpansion2CK(state::Dict{BitVector,Float64}, sectors::String)
    @assert sectors in ["p", "h", "ph"] "Provided IOM sectors not among p, h or ph/hp."
    if sectors == "p"
        state = Dict{keytype(state),valtype(state)}([key; [1, 1, 1, 1]] => val for (key, val) in state)
    elseif sectors == "h"
        state = Dict{keytype(state),valtype(state)}([key; [0, 0, 0, 0]] => val for (key, val) in state)
    else
        state = Dict{keytype(state),valtype(state)}([key; [1, 1, 0, 0, 1, 1, 0, 0]] => val for (key, val) in state)
    end
end


function stateExpansion1CK(state::Dict{BitVector,Float64}, sectors::String)
    @assert sectors in ["p", "h", "ph"] "Provided IOM sectors not among p, h or ph/hp."
    if sectors == "p"
        state = Dict{keytype(state),valtype(state)}([key; [1, 1]] => val for (key, val) in state)
    elseif sectors == "h"
        state = Dict{keytype(state),valtype(state)}([key; [0, 0]] => val for (key, val) in state)
    else
        state = Dict{keytype(state),valtype(state)}([key; [1, 1, 0, 0]] => val for (key, val) in state)
    end
    return state
end
