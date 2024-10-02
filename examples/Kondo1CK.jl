using fermions, Plots, Measures
theme(:dark)
include("../src/iterDiag.jl")

totalSites = 8
initSites = 1
kondoJ = 1.
maxSize = 500
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
p1 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="\$\\langle S_d^+ S_{k_1}^- + \\mathrm{h.c.} \\rangle\$")
p2 = plot(thickness_scaling=1.6, leftmargin=-6mm, bottommargin=-3mm, label="approx.", xlabel="sites", ylabel="energy per site")
spinFlipCorrd2 = Tuple{String, Vector{Int64}, Float64}[("+-+-", [1, 2, 8, 7], 1.0), ("+-+-", [2, 1, 7, 8], 1.0)]
spinFlipCorrd0 = Tuple{String, Vector{Int64}, Float64}[("+-+-", [1, 2, 4, 3], 1.0), ("+-+-", [2, 1, 3, 4], 1.0)]
savePaths, resultsDict = IterDiag(hamFlow, maxSize;
                     symmetries=Char['N', 'S'],
                     # symmetries=Char['S'],
                     # symmetries=Char['N'],
                     # magzReq=(m, N) -> -2 ≤ m ≤ 2,
                     # occReq=(x, N) -> div(N, 2) - 4 ≤ x ≤ div(N, 2) + 4,
                     correlationDefDict=Dict("SF0" => spinFlipCorrd0, "SF2" => spinFlipCorrd2),
                    ) 
# display(deserialize(savePaths[end-1])["eigVals"])
# display(deserialize(savePaths[end-1])["results"]["SF2"])
scatter!(p1, [deserialize(path)["results"]["SF2"] for path in savePaths[1:end-1] if !isnothing(deserialize(path)["results"]["SF2"])])
scatter!(p2, initSites:totalSites, [deserialize(path)["eigVals"][1] / (initSites + i) for (i, path) in enumerate(savePaths[1:end-1])])

corrExact = []
energyExact = []
@showprogress for (i, num) in enumerate(initSites:totalSites)
    basis = BasisStates(2 * (1 + num); totOccReq=[1 + num], magzReq=[ifelse(isodd(i), 0, 1)], localCriteria=x->x[1]+x[2]==1)
    fullHam = vcat(hamFlow[1:i]...)
    fullMatrix = OperatorMatrix(basis, fullHam)
    F = eigen(fullMatrix)
    push!(energyExact, minimum(F.values)/(1 + num))
    try
        totalNum = OperatorMatrix(basis, [("n", [i], 1.) for i in 1:2*(1+num)])
        push!(corrExact, F.vectors[:, 1]' * OperatorMatrix(basis, spinFlipCorrd2) * F.vectors[:, 1])
    catch e
        continue
    end
end
display(energyExact)
display(corrExact)
# 
plot!(p1, corrExact, label="exact", linewidth=2)
plot!(p2, initSites:totalSites, energyExact, label="exact", linewidth=2)
savefig(p1, "spinflipcomparison.pdf")
savefig(p2, "energy.pdf")
