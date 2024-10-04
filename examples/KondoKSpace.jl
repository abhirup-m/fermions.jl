using fermions, Plots, Measures, BenchmarkTools
theme(:dark)
include("../src/iterDiag.jl")

totalSites = 7
initSites = 1
kondoJ = 2.
bathInt = -1.5
maxSize = 2900
dispersion = [ifelse(i % 2 == 1, -1., 1.) for i in 1:totalSites]

spinFlipCorrd2 = Tuple{String, Vector{Int64}, Float64}[("+-+-", [1, 2, 2 * totalSites + 2, 2 * totalSites + 1], 1.0), ("+-+-", [2, 1, 2 * totalSites + 1, 2 * totalSites + 2], 1.0)]
spinFlipCorrd0 = Tuple{String, Vector{Int64}, Float64}[("+-+-", [1, 2, 4, 3], 1.0), ("+-+-", [2, 1, 3, 4], 1.0)]
vneDefDict = Dict("vne_d" => [1, 2], "vne_0" => [5, 6], "vne_d0" => [1, 2, 7, 8])
mutInfoDefDict = Dict("I2_d0" => ([2 * totalSites + 1, 2 * totalSites + 2], [3, 4]))
hamFlow = getHamFlow(initSites, totalSites, hop_t, kondoJ)

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

    for site in initSites+1:totalSites
        newTerm = []
        push!(newTerm, ("n",  [1 + 2 * site], dispersion[site])) # c^†_{j,up} c_{j+1,up}
        push!(newTerm, ("n",  [2 + 2 * site], dispersion[site])) # c^†_{j+1,up} c_{j,up}

        for (k1, k2) in Iterators.product(1:site, 1:site)
            if k1 ≠ site && k2 ≠ site
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
    p1 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="\$\\langle S_d^+ S_{k_1}^- + \\mathrm{h.c.} \\rangle\$")
    p2 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="energy per site")
    savePaths, resultsDict = IterDiag(hamFlow, maxSize;
                         symmetries=Char['N', 'S'],
                         # symmetries=Char['S'],
                         #=symmetries=Char['N'],=#
                         magzReq=(m, N) -> -2 ≤ m ≤ 3,
                         occReq=(x, N) -> div(N, 2) - 3 ≤ x ≤ div(N, 2) + 3,
                         corrMagzReq=(m, N) -> m == ifelse(isodd(div(N, 2)), 1, 0),
                         corrOccReq=(x, N) -> x == div(N, 2),
                         correlationDefDict=Dict("SF0" => spinFlipCorrd0, "SF2" => spinFlipCorrd2),
                         vneDefDict=deepcopy(vneDefDict),
                         mutInfoDefDict=deepcopy(mutInfoDefDict),
                        ) 
    display(resultsDict)
    #=@assert false=#
end

function ExactResults()
    corrExact = []
    energyExact = []
    vneExact = []
    mutInfoExact = []
    finalState = nothing
    finalEnergy = 0
    @showprogress for (i, num) in enumerate(initSites:totalSites)
        basis = BasisStates(2 * (1 + num); totOccReq=[1 + num], magzReq=[ifelse(isodd(i), 0, 1)], localCriteria=x->x[1]+x[2]==1)
        fullHam = vcat(hamFlow[1:i]...)
        E, X = Spectrum(fullHam, basis)
        finalState = X[1]
        finalEnergy = E[1] / (2 * (1 + num))
    end
    display(finalEnergy)
    display(GenCorrelation(finalState, spinFlipCorrd2))
    for (name, sites) in vneDefDict
        println(name, ": ", VonNEntropy(finalState, sites))
    end
    for (name, sites) in mutInfoDefDict
        println(name, ": ", MutInfo(finalState, sites))
    end
    # 
    #=plot!(p1, corrExact, label="exact", linewidth=2)=#
    #=plot!(p2, initSites:totalSites, energyExact, label="exact", linewidth=2)=#
    #=savefig(p1, "spinflipcomparison.pdf")=#
    #=savefig(p2, "energy.pdf")=#
end
IterResults()
ExactResults()
