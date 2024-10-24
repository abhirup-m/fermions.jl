using fermions, Plots, Measures
include("../src/iterDiag.jl")
include("../src/correlations.jl")

totalSites = 2
initSites = 1
maxSize = 5000
kondoJ = 1.
hop_t = 1.1

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

    for site in initSites+1:totalSites
        newTerm = []
        push!(newTerm, ("+-",  [1 + 2 * site, -1 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
        push!(newTerm, ("+-",  [-1 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
        push!(newTerm, ("+-",  [2 + 2 * site, 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
        push!(newTerm, ("+-",  [2 * site, 2 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}

        push!(hamFlow, newTerm)
    end
    return hamFlow
end

hamFlow = getHamFlow(initSites, totalSites, hop_t, kondoJ)

savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                             symmetries=Char['N', 'S'],
                             # symmetries=Char['S'],
                             #=symmetries=Char['N'],=#
                             #=magzReq=(m, N) -> -2 ≤ m ≤ 3,=#
                             #=occReq=(x, N) -> div(N, 2) - 4 ≤ x ≤ div(N, 2) + 4,=#
                             #=corrMagzReq=(m, N) -> m == ifelse(isodd(div(N, 2)), 1, 0),=#
                             #=corrOccReq=(x, N) -> x == div(N, 2),=#
                             specFuncDefDict=Dict("create" => ("+", [1]), "destroy" => ("-", [1])),
                            )
#=specFuncOperators = Dict("destroy" => specFuncOperators["imp"][1], "create" => specFuncOperators["imp"][2])=#
freqValues = collect(-1:0.5:1)
standDev = 0.02
@time specFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                        occReq=(x,N) -> x == div(N,2),
                        #=excOccReq=(x,N) -> abs(x - div(N,2)) == 1,=#
                       )
#=p = plot(freqValues, specFunc)=#
#=savefig(p, "specfunc.pdf")=#
#=display(p)=#


function ExactResults()
    totalSpecFunc = zeros(length(freqValues))
    @showprogress for (i, num) in enumerate(initSites:totalSites)
        basis = BasisStates(2 * (1 + num))#; totOccReq=[num, 1 + num, 2 + num], magzReq=[-1, 0, 1])#, localCriteria=x->x[1]+x[2]==1)
        fullHam = vcat(hamFlow[1:i]...)
        E, X = Spectrum(fullHam, basis)
        specFunc = SpecFunc(E, X, [("-", [1], 1.)], [("+", [1], 1.)], freqValues, basis, standDev)
        println(round.(specFunc, digits=5))
        totalSpecFunc .+= specFunc
    end
    println(round.(totalSpecFunc, digits=5))
end
ExactResults()
