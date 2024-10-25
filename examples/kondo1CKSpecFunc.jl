using fermions, Plots, Measures
include("../src/iterDiag.jl")

initSites = 1
maxSize = 1500
kondoJ = 1.
hop_t = 1.

function getHamFlow(initSites::Int64, totalSites::Int64, hop_t::Float64, kondoJ::Float64)
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[]
    initHam = Tuple{String, Vector{Int64}, Float64}[]
    for site in 1:(initSites-1)
        push!(initHam, ("+-",  [1 + 2 * site, 3 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
        push!(initHam, ("+-",  [3 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
        push!(initHam, ("+-",  [2 + 2 * site, 4 + 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
        push!(initHam, ("+-",  [4 + 2 * site, 2 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}
    end
    push!(initHam, ("nn",  [1, 3], kondoJ/4)) # n_{d up, n_{0 up}
    push!(initHam, ("nn",  [1, 4], -kondoJ/4)) # n_{d up, n_{0 down}
    push!(initHam, ("nn",  [2, 3], -kondoJ/4)) # n_{d down, n_{0 up}
    push!(initHam, ("nn",  [2, 4], kondoJ/4)) # n_{d down, n_{0 down}
    push!(initHam, ("+-+-",  [1, 2, 4, 3], kondoJ/2)) # S_d^+ S_0^-
    push!(initHam, ("+-+-",  [2, 1, 3, 4], kondoJ/2)) # S_d^- S_0^+

    push!(hamFlow, initHam)

    for site in initSites+1:2:totalSites
        newTerm = []
        push!(newTerm, ("+-",  [1 + 2 * site, -1 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
        push!(newTerm, ("+-",  [-1 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
        push!(newTerm, ("+-",  [2 + 2 * site, 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
        push!(newTerm, ("+-",  [2 * site, 2 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}

        push!(newTerm, ("+-",  [3 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
        push!(newTerm, ("+-",  [1 + 2 * site, 3 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
        push!(newTerm, ("+-",  [4 + 2 * site, 2 + 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
        push!(newTerm, ("+-",  [2 + 2 * site, 4 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}

        push!(hamFlow, newTerm)
    end
    return hamFlow
end

p = plot()
for totalSites in 20:20:100
    hamFlow = getHamFlow(initSites, totalSites, hop_t, kondoJ)
    specFuncDefDict = Dict("create" => ("+-+", [1,2,4]), "destroy" => ("+--", [2,1,4]))
    savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                                 symmetries=Char['N', 'S'],
                                 specFuncDefDict=specFuncDefDict,
                                 occReq=(x,N)->abs(x-div(N,2)) ≤ 4
                                )
    freqValues = collect(-4:0.01:4)
    standDev = 0.1
    @time totalSpecFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                                       occReq=(x,N) -> x == div(N,2),
                           )
    plot!(p, freqValues, totalSpecFunc)
end


#=function ExactResults()=#
#=    totalSpecFunc = zeros(length(freqValues))=#
#=    @showprogress for (i, num) in enumerate(initSites:2:totalSites)=#
#=        basis = BasisStates(2 * (1 + num); localCriteria=x->x[1]+x[2]==1, totOccReq=[num, 1 + num, 2 + num])=#
#=        fullHam = vcat(hamFlow[1:i]...)=#
#=        E, X = Spectrum(fullHam, basis)=#
#=        specFunc = SpecFunc(E, X, Dict("create" => [(specFuncDefDict["create"]..., 1.)], "destroy" => [(specFuncDefDict["destroy"]..., 1.)]), =#
#=                            freqValues, basis, standDev, ['N'], (1+num,))=#
#=        totalSpecFunc .+= specFunc=#
#=    end=#
#=    return totalSpecFunc=#
#=end=#
#=totalSpecFunc = ExactResults()=#
#=plot!(p, freqValues, totalSpecFunc, linestyle=:dot)=#
display(p)
