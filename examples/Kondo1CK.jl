using fermions, Plots, Measures, BenchmarkTools
theme(:dark)
include("../src/iterDiag.jl")

totalSites = 12
initSites = 1
kondoJ = abs(rand())
bathInt = -1.5
maxSize = 200
dispersion = sort(rand(totalSites) .^ 2)

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

    # push!(initHam, ("n",  [3], -bathInt/2)) # -Ub/2 n_{0 up}
    # push!(initHam, ("n",  [4], -bathInt/2)) # -Ub/2 n_{0 down}
    # push!(initHam, ("nn",  [3, 4], bathInt)) # Ub n_{0 up} n_{0 down}

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
            if k2 == site
                push!(newTerm, ("n+-",  [1, up1, up2], kondoJ/4)) # n_{d up, n_{0 up}
                push!(newTerm, ("n+-",  [1, down1, down2], -kondoJ/4)) # n_{d up, n_{0 down}
                push!(newTerm, ("n+-",  [2, up1, up2], -kondoJ/4)) # n_{d down, n_{0 up}
                push!(newTerm, ("n+-",  [2, down1, down2], kondoJ/4)) # n_{d down, n_{0 down}
                push!(newTerm, ("+-+-",  [1, 2, down1, up2], kondoJ/2)) # S_d^+ S_0^-
                push!(newTerm, ("+-+-",  [2, 1, up1, down2], kondoJ/2)) # S_d^- S_0^+
            else
                push!(newTerm, ("+-n",  [up1, up2, 1], kondoJ/4)) # n_{d up, n_{0 up}
                push!(newTerm, ("+-n",  [down1, down2, 1], -kondoJ/4)) # n_{d up, n_{0 down}
                push!(newTerm, ("+-n",  [up1, up2, 2], -kondoJ/4)) # n_{d down, n_{0 up}
                push!(newTerm, ("+-n",  [down1, down2, 2], kondoJ/4)) # n_{d down, n_{0 down}
                push!(newTerm, ("+-+-",  [down1, up2, 1, 2], kondoJ/2)) # S_d^+ S_0^-
                push!(newTerm, ("+-+-",  [up1, down2, 2, 1], kondoJ/2)) # S_d^- S_0^+
            end
        end

        push!(hamFlow, newTerm)
    end
    return hamFlow
end

hamFlow = getHamFlow(initSites, totalSites, dispersion, kondoJ, bathInt)
p1 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="\$\\langle S_d^+ S_{k_1}^- + \\mathrm{h.c.} \\rangle\$")
p2 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="energy per site")
corrdef = [("+-+-", [1, 2, 4, 3], 1.0), ("+-+-", [2, 1, 3, 4], 1.0)]
savePaths = IterDiag(hamFlow, maxSize; symmetries=Char['N'], occReq=(x, N) -> ifelse(0 ≤ x ≤ 4, true, false))
corrFlow, energypersite = IterCorrelation(savePaths, corrdef; occupancy=(x, N) -> x == 2)
display(energypersite)
display(corrFlow)
scatter!(p1, initSites:totalSites, corrFlow, label="\$M_s=$(maxSize)\$")
scatter!(p2, initSites:totalSites, energypersite, label="\$M_s=$(maxSize)\$")

corrExact = []
energyExact = []
@showprogress for (i, num) in enumerate(initSites:totalSites)
    basis = BasisStates(2 * (1 + num); totOccReq=[2], magzReq=[0], localCriteria=x->x[1]+x[2]==1)
    fullHam = vcat(hamFlow[1:i]...)
    fullMatrix = OperatorMatrix(basis, fullHam)
    # display(OperatorMatrix(basis, [("+h", [2, 1], 1.)]))
    # println("----")
    F = eigen(fullMatrix)
    push!(corrExact, F.vectors[:, 1]' * OperatorMatrix(basis, corrdef) * F.vectors[:, 1])
    push!(energyExact, minimum(F.values)/(2*(1 + num)))
end
display(energyExact)
display(corrExact)
# 
plot!(p1, initSites:totalSites, corrExact, label="exact", linewidth=2)
plot!(p2, initSites:totalSites, energyExact, label="exact", linewidth=2)
savefig(p1, "spinflipcomparison.pdf")
savefig(p2, "energy.pdf")
