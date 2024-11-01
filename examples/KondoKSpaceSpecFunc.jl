using fermions, Plots, Measures, BenchmarkTools
theme(:dark)
include("../src/iterDiag.jl")

initSites = 1
kondoJ = 2.
bathInt = -1.5
maxSize = 2000

#=spinFlipCorrd2 = Tuple{String, Vector{Int64}, Float64}[("+-+-", [1, 2, 2 * totalSites + 2, 2 * totalSites + 1], 1.0), ("+-+-", [2, 1, 2 * totalSites + 1, 2 * totalSites + 2], 1.0)]=#
#=spinFlipCorrd0 = Tuple{String, Vector{Int64}, Float64}[("+-+-", [1, 2, 4, 3], 1.0), ("+-+-", [2, 1, 3, 4], 1.0)]=#
#=vneDefDict = Dict("vne_d" => [1, 2], "vne_0" => [5, 6], "vne_d0" => [1, 2, 7, 8])=#
#=mutInfoDefDict = Dict("I2_d0" => ([2 * totalSites + 1, 2 * totalSites + 2], [3, 4]))=#
#=hamFlow = getHamFlow(initSites, totalSites, hop_t, kondoJ)=#

function getHamFlow(initSites::Int64, totalSites::Int64, dispersion::Vector{Float64}, kondoJ::Float64, bathInt::Float64)
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[]
    initHam = Tuple{String, Vector{Int64}, Float64}[]
    for site in 1:initSites
        push!(initHam, ("n",  [1 + 2 * site], dispersion[site])) # c^†_{j,up} c_{j+1,up}
        push!(initHam, ("n",  [2 + 2 * site], dispersion[site])) # c^†_{j+1,up} c_{j,up}
    end
    for (k1, k2) in Iterators.product(1:initSites, 1:initSites)
        up1, up2 = 2 .* [k1, k2] .+ 1
        down1, down2 = 2 .* [k1, k2] .+ 2
        push!(initHam, ("n+-",  [1, up1, up2], kondoJ/4)) # n_{d up, n_{0 up}
        push!(initHam, ("n+-",  [1, down1, down2], -kondoJ/4)) # n_{d up, n_{0 down}
        push!(initHam, ("n+-",  [2, up1, up2], -kondoJ/4)) # n_{d down, n_{0 up}
        push!(initHam, ("n+-",  [2, down1, down2], kondoJ/4)) # n_{d down, n_{0 down}
        push!(initHam, ("+-+-",  [1, 2, down1, up2], kondoJ/2)) # S_d^+ S_0^-
        push!(initHam, ("+-+-",  [2, 1, up1, down2], kondoJ/2)) # S_d^- S_0^+
    end

    push!(hamFlow, initHam)

    for site in initSites+1:2:totalSites
        newTerm = []
        push!(newTerm, ("n",  [1 + 2 * site], dispersion[site])) # c^†_{j,up} c_{j+1,up}
        push!(newTerm, ("n",  [2 + 2 * site], dispersion[site])) # c^†_{j+1,up} c_{j,up}
        push!(newTerm, ("n",  [3 + 2 * site], dispersion[site+1])) # c^†_{j,up} c_{j+1,up}
        push!(newTerm, ("n",  [4 + 2 * site], dispersion[site+1])) # c^†_{j+1,up} c_{j,up}

        for (k1, k2) in Iterators.product(1:site+1, 1:site+1)
            if k1 ∉ (site, site+1) && k2 ∉ (site, site+1)
                continue
            end
            up1, up2 = 2 .* [k1, k2] .+ 1
            down1, down2 = 2 .* [k1, k2] .+ 2
            push!(newTerm, ("n+-",  [1, up1, up2], kondoJ/4)) # n_{d up, n_{0 up}
            push!(newTerm, ("n+-",  [1, down1, down2], -kondoJ/4)) # n_{d up, n_{0 down}
            push!(newTerm, ("n+-",  [2, up1, up2], -kondoJ/4)) # n_{d down, n_{0 up}
            push!(newTerm, ("n+-",  [2, down1, down2], kondoJ/4)) # n_{d down, n_{0 down}
            push!(newTerm, ("+-+-",  [1, 2, down1, up2], kondoJ/2)) # S_d^+ S_0^-
            push!(newTerm, ("+-+-",  [2, 1, up1, down2], kondoJ/2)) # S_d^- S_0^+
        end

        push!(hamFlow, newTerm)
    end
    return hamFlow
end

function IterResults()
    p = plot()
    for totalSites in 3:2:7
        dispersion = [0. for i in 1:totalSites]
        hamFlow = getHamFlow(initSites, totalSites, dispersion, kondoJ, bathInt)
        specFuncDefDict = Dict("create" => [("+-+", [1, 2, 4], 1.)], "destroy" => [("+--", [2, 1, 4], 1.)])
        #=specFuncDefDict = Dict("create" => [("+-+", [1, 2, 2 + 2 * i], 1.) for i in 1:totalSites], "destroy" => [("+--", [2, 1, 2 + 2 * i], 1.) for i in 1:totalSites])=#
        savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                                     symmetries=Char['N', 'S'],
                                     specFuncDefDict=specFuncDefDict,
                                     occReq=(x,N)->abs(x-div(N,2)) ≤ 4
                                    )
        freqValues = collect(-6:0.01:6)
        standDev = 0.1
        @time totalSpecFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                                           occReq=(x,N) -> x == div(N,2),
                               )
        plot!(p, freqValues, totalSpecFunc)
    end
    display(p)
end

function ExactResults()
    freqValues = collect(-6:0.01:6)
    standDev = 0.1
    p = plot()
    for totalSites in 3:2:5
        specFuncDefDict = Dict("create" => [("+-+", [1, 2, 2 + 2 * i], 1.) for i in 1:totalSites], "destroy" => [("+--", [2, 1, 2 + 2 * i], 1.) for i in 1:totalSites])
        dispersion = [0. for i in 1:totalSites]
        totalSpecFunc = zeros(length(freqValues))
        @showprogress for (i, num) in enumerate(initSites:2:totalSites)
            basis = BasisStates(2 * (1 + num); localCriteria=x->x[1]+x[2]==1, totOccReq=[num, 1 + num, 2 + num])
            hamFlow = getHamFlow(initSites, totalSites, dispersion, kondoJ, bathInt)
            fullHam = vcat(hamFlow[1:i]...)
            println(hamFlow[end])
            E, X = Spectrum(fullHam, basis)
            specFunc = SpecFunc(E, X, Dict("create" => specFuncDefDict["create"], "destroy" => specFuncDefDict["destroy"]), 
                                freqValues, basis, standDev, ['N'], (1+num,))
            totalSpecFunc .+= specFunc
        end
        plot!(p, freqValues, totalSpecFunc, linestyle=:dot)
    end
    display(p)
end
IterResults()
#=ExactResults()=#
