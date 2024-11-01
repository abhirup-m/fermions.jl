using fermions, CairoMakie, Measures, ProgressMeter
include("../src/base.jl")
include("../src/correlations.jl")
include("../src/iterDiag.jl")

set_theme!(merge(theme_light(), theme_latexfonts()))
update_theme!(
              figure_padding = 0,
              fontsize=28,
              ScatterLines = (
                       linewidth = 3,
                       markersize=10,
                      ),
              Lines = (
                       linewidth = 3,
                       markersize=20,
                      ),
              Scatter = (
                       markersize=10,
                      ),
              Legend = (
                        patchsize=(50,20),
                        halign = :right,
                        valign = :top,
                       ),
             )

function getHamFlow(initSites::Int64, totalSites::Int64, hop_t::Float64, hyb::Float64, impCorr::Float64, field::Float64,)
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[]
    initHam = Tuple{String, Vector{Int64}, Float64}[]
    for site in 1:(initSites-1)
        push!(initHam, ("+-",  [1 + 2 * site, 3 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
        push!(initHam, ("+-",  [3 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
        push!(initHam, ("+-",  [2 + 2 * site, 4 + 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
        push!(initHam, ("+-",  [4 + 2 * site, 2 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}
    end

    for site in 0:initSites
        push!(initHam, ("n",  [1 + 2 * site], field))
        push!(initHam, ("n",  [2 + 2 * site], -field))
    end
    push!(initHam, ("n",  [1], -impCorr/2))
    push!(initHam, ("n",  [2], -impCorr/2))
    push!(initHam, ("nn",  [1,2], impCorr))

    push!(initHam, ("+-",  [1, 3], -hyb)) # c^†_{j,up} c_{j+1,up}
    push!(initHam, ("+-",  [3, 1], -hyb)) # c^†_{j+1,up} c_{j,up}
    push!(initHam, ("+-",  [2, 4], -hyb)) # c^†_{j,dn} c_{j+1,dn}
    push!(initHam, ("+-",  [4, 2], -hyb)) # c^†_{j+1,dn} c_{j,dn}

    push!(hamFlow, initHam)

    for site in initSites+1:addPerStep:totalSites
        newTerm = []
        for newSite in site:(site+addPerStep-1)
            push!(newTerm, ("+-",  [1 + 2 * newSite, -1 + 2 * newSite], -hop_t)) # c^†_{j,up} c_{j+1,up}
            push!(newTerm, ("+-",  [-1 + 2 * newSite, 1 + 2 * newSite], -hop_t)) # c^†_{j+1,up} c_{j,up}
            push!(newTerm, ("+-",  [2 + 2 * newSite, 2 * newSite], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
            push!(newTerm, ("+-",  [2 * newSite, 2 + 2 * newSite], -hop_t)) # c^†_{j+1,dn} c_{j,dn}

            push!(newTerm, ("n",  [1 + 2 * newSite], field))
            push!(newTerm, ("n",  [-1 + 2 * newSite], -field))
        end

        push!(hamFlow, newTerm)
    end
    return hamFlow
end

function IterResults(hamFlow, totalSites::Int64)
    savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                                                         symmetries=Char['N'],#, 'S'],
                                 #=correlationDefDict=Dict("ndup" => [("n", [1], 1.)], "nddown" => [("n", [2], 1.)]),=#
                                 specFuncDefDict=specFuncDefDict,
                                 #=occReq=(x,N)->abs(x-div(N,2)) ≤ 4=#
                                )
    totalSpecFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                                       occReq=(x,N) -> x == div(N,2),
                                       excOccReq=(x,N) -> abs(x - div(N,2)) == 1,
                           )
    return totalSpecFunc
end

function ExactResults(hamFlow, totalSites::Int64)
    totalSpecFunc = zeros(length(freqValues))
    for (i, num) in enumerate(initSites:addPerStep:totalSites)
        basis = BasisStates(2 * (1 + num); 
                            #=localCriteria=x->x[1]+x[2]==1,=#
                            totOccReq=[num, 1 + num, 2 + num]
                           )
        fullHam = vcat(hamFlow[1:i]...)
        E, X = Spectrum(fullHam, basis)
        specFunc = SpecFunc(E, X,
                            Dict("create" => specFuncDefDict["create"],
                                 "destroy" => specFuncDefDict["destroy"]), 
                            freqValues, basis, standDev, 
                            ['N'], (1+num,)
                           )
        totalSpecFunc .+= specFunc
    end
    return totalSpecFunc
end

function BenchMark(hop_t, hyb, impCorr)
    f = Figure()
    ax = Axis(f[1,1], xlabel=L"\omega",ylabel=L"A(\omega)")
    for totalSites in [5]
        hamFlow = getHamFlow(initSites, totalSites, hop_t, hyb, impCorr, 1e-5)
        specFuncIter = IterResults(hamFlow, totalSites)
        #=println(round.(specFuncIter, digits=10))=#
        #=println("------------------")=#
        hamFlow = getHamFlow(initSites, totalSites, hop_t, hyb, impCorr, 1e-5)
        specFuncExact = ExactResults(hamFlow, totalSites)
        #=println(round.(specFuncIter, digits=10))=#
        lines!(ax, freqValues, specFuncIter, label=L"ID $L=%$(totalSites)$")
        scatter!(ax, freqValues, specFuncExact, linestyle=:dash, label=L"ED $L=%$(totalSites)$")
    end
    axislegend(ax)
    #=display(f)=#
    save("specFuncComparison.pdf", f)
end


function LargerSystem(hop_t, hyb, impCorr, standDev)
    f = Figure()
    ax = Axis(f[1,1], xlabel=L"\omega",ylabel=L"A(\omega)", title=L"\eta=%$(maximum(standDev)),t=%$(hop_t)")#, yscale=log10)
    totalSites = 21
    for hyb in [0.1, 0.0001, 10^(-5.2), 10^(-5.3), 0.]
        hamFlow = getHamFlow(initSites, totalSites, hop_t, hyb, impCorr, -1e-5)
        specFuncIter = IterResults(hamFlow, totalSites)
        specFuncIter ./= sum(specFuncIter .* (maximum(freqValues) - minimum(freqValues[1])) / (length(freqValues) - 1))
        scatterlines!(ax, freqValues, 1e-5 .+ specFuncIter, label=L"$V/U=%$(abs(hyb)/impCorr)$")
    end
    axislegend(ax)
    #=display(f)=#
    save("specFunc-$(maximum(standDev))-gap.pdf", f)
    #=save("specFunc-$(standDev)-$(hop_t).pdf", f)=#
end

specFuncDefDict = Dict("create" => [("+", [2], 1.), ("+", [1], 1.)], "destroy" => [("-", [2], 1.), ("-", [1], 1.)])
initSites = 1
maxSize = 1000
hop_t = 0.1
hyb = 0.1
impCorr = 4.
addPerStep = 1
standDev = ones(length(freqValues))
standDev[freqValues .>= 0] .= 0.08 # range(0.01, stop=0.1, length=sum(freqValues .>= 0))
standDev[freqValues .<= 0] .= 0.08 # range(0.1, stop=0.01, length=sum(freqValues .<= 0))
println((maximum(standDev), minimum(standDev)))
freqValues = collect(-3:0.01:3)
LargerSystem(hop_t, hyb, impCorr, standDev)
#=BenchMark(kondoJ, hop_t)=#
