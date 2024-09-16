using fermions, Plots, Measures, BenchmarkTools
theme(:dark)
include("../src/iterDiag.jl")

realSites = 2
initSites = 1
dispersion = -1. .* collect(1:realSites)
kondoJMatrix = 2.1 .* ones(realSites, realSites)
bathIntMatrix = 0.0 * ones(realSites, realSites, realSites, realSites)

function getHamFlow(initSites::Int64, realSites::Int64, dispersion::Vector{Float64}, kondoJMatrix::Matrix{Float64}, bathIntMatrix::Array{Float64, 4})
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[]
    initHam = []
    for k in 1:initSites
        up = 1 + 2 * k
        down = up + 1
        push!(initHam, ("n",  [up], dispersion[k]))
        push!(initHam, ("n",  [down], dispersion[k]))
    end
    for (k1, k2) in Iterators.product(1:initSites, 1:initSites)
        up1 = 1 + 2 * k1
        up2 = 1 + 2 * k2
        down1, down2 = up1 + 1, up2 + 1

        kondoJ = kondoJMatrix[k1, k2]
        # push!(initHam, ("n+-",  [1, up1, up2], kondoJ/4))
        # push!(initHam, ("n+-",  [1, down1, down2], -kondoJ/4))
        # push!(initHam, ("n+-",  [2, up1, up2], -kondoJ/4))
        # push!(initHam, ("n+-",  [2, down1, down2], kondoJ/4))
        push!(initHam, ("+-+-",  [1, 2, down1, up2], kondoJ/2))
        push!(initHam, ("+-+-",  [2, 1, up1, down2], kondoJ/2))
    end
    # for k_i in Iterators.product(repeat([1:initSites], 4)...)
    #     up1, up2, up3, up4 = 1 .+ 2 .* k_i
    #     down1, down2, down3, down4 = [up1, up2, up3, up4] .+ 1
    #     bathInt = bathIntMatrix[k_i...]
    #     push!(initHam, ("+-+-",  [up1, up2, up3, up4], -bathInt/2))
    #     push!(initHam, ("+-+-",  [down1, down2, down3, down4], -bathInt/2))
    #     push!(initHam, ("+-+-",  [up1, up2, down3, down4], bathInt/2))
    #     push!(initHam, ("+-+-",  [down1, down2, up3, up4], bathInt/2))
    # end

    push!(hamFlow, initHam)
    for site in initSites+1:realSites
        newTerm = []
        push!(newTerm, ("n", [1 + 2 * site], dispersion[site]))
        push!(newTerm, ("n", [2 + 2 * site], dispersion[site]))
        # push!(newTerm, ("n", [3 + 2 * site], dispersion[site+1]))
        # push!(newTerm, ("n", [4 + 2 * site], dispersion[site+1]))
        for (k1, k2) in Iterators.product(1:site, 1:site)
            if site ∉ (k1, k2)# || site + 1 ∈ (k1, k2)
                continue
            end
            up1 = 1 + 2 * k1
            up2 = 1 + 2 * k2
            kondoJ = kondoJMatrix[k1, k2]
            down1, down2 = up1 + 1, up2 + 1
            # push!(newTerm, ("n+-",  [1, up1, up2], kondoJ/4))
            # push!(newTerm, ("n+-",  [1, down1, down2], -kondoJ/4))
            # push!(newTerm, ("n+-",  [2, up1, up2], -kondoJ/4))
            # push!(newTerm, ("n+-",  [2, down1, down2], kondoJ/4))
            push!(newTerm, ("+-+-",  [1, 2, down1, up2], kondoJ/2))
            push!(newTerm, ("+-+-",  [2, 1, up1, down2], kondoJ/2))
        end
        # for k_i in Iterators.product(repeat([1:site], 4)...)
        #     if site ∉ k_i || site + 1 ∈ k_i
        #         continue
        #     end
        #     up1, up2, up3, up4 = 1 .+ 2 .* k_i
        #     down1, down2, down3, down4 = [up1, up2, up3, up4] .+ 1
        #     bathInt = bathIntMatrix[k_i...]
        #     push!(initHam, ("+-+-",  [up1, up2, up3, up4], -bathInt/2))
        #     push!(initHam, ("+-+-",  [down1, down2, down3, down4], -bathInt/2))
        #     push!(initHam, ("+-+-",  [up1, up2, down3, down4], bathInt/2))
        #     push!(initHam, ("+-+-",  [down1, down2, up3, up4], bathInt/2))
        # end
        push!(hamFlow, newTerm)
    end
    return hamFlow
end

hamFlow = getHamFlow(initSites, realSites, dispersion, kondoJMatrix, bathIntMatrix)
p1 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="\$\\langle S_d^+ S_{k_1}^- + \\mathrm{h.c.} \\rangle\$")
p2 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="energy per site")
corrdef = [("+-+-", [1, 2, 4, 3], 1.0), ("+-+-", [2, 1, 3, 4], 1.0)]
for maxSize in [600]
    savePaths = IterDiag(hamFlow, maxSize; symmetries=Char['N'])#, occReq=(x, N) -> ifelse(N-4 ≤ x ≤ N+4, true, false))
    display(minimum(deserialize(savePaths[end-1])["eigVals"][deserialize(savePaths[end-1])["occupancy"] .== 3]))
    corrFlow, energypersite = IterCorrelation(savePaths, corrdef)
    # scatter!(p1, initSites:2:realSites, corrFlow, label="\$M_s=$(maxSize)\$")
    # scatter!(p2, initSites:2:realSites, energypersite .* (1 .+ (initSites:2:realSites)), label="\$M_s=$(maxSize)\$")
end

corrExact = []
energyExact = []
@showprogress for (i, num) in enumerate(initSites:realSites)
    basis = BasisStates(2 * (1 + num); totOccReq=[1 + num])
    fullHam = vcat(hamFlow[1:i]...)
    fullMatrix = OperatorMatrix(basis, fullHam)
    display(fullMatrix)
    F = eigen(fullMatrix)
    # push!(corrExact, F.vectors[:, 1]' * OperatorMatrix(basis, corrdef) * F.vectors[:, 1])
    # push!(energyExact, minimum(F.values))
    println(F.values[1])
end

# plot!(p1, initSites:2:realSites, corrExact, label="exact", linewidth=2)
# plot!(p2, initSites:2:realSites, energyExact, label="exact", linewidth=2)
# savefig(p1, "spinflipcomparison.pdf")
# savefig(p2, "energy.pdf")
