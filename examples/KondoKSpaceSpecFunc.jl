using fermions, CairoMakie, Measures, ProgressMeter
include("../src/iterDiag.jl")
include("../src/correlations.jl")

set_theme!(merge(theme_ggplot2(), theme_latexfonts()))
update_theme!(
              figure_padding = 0,
              fontsize=28,
              ScatterLines = (
                       linewidth = 3,
                       markersize=10,
                      ),
              Lines = (
                       linewidth = 6,
                       markersize=20,
                      ),
              Scatter = (
                       markersize=8,
                      ),
              Legend = (
                        patchsize=(50,20),
                        halign = :right,
                        valign = :top,
                       ),
             )


function IterResults(hamFlow, totalSites::Int64, specFuncDefDict, maxSize::Int64)
    savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                                                         symmetries=Char['N'],#, 'S'],
                                                         specFuncDefDict=specFuncDefDict,
                                                         #=occReq=(x,N)->abs(x-div(N,2)) ≤ 4=#
                                                         )
    totalSpecFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                                 degenTol=1e-10,
                                 occReq=(x,N) -> x == div(N,2),
                                 excOccReq=(x,N) -> abs(x - div(N,2)) == 1,
                                 #=symmetrise=true,=#
                           )
    return totalSpecFunc
end


function LargerSystem(
        totalSites::Int64,
        hamFlow;
    )
    f = Figure()
    ax = Axis(f[1,1], xlabel=L"\omega",ylabel=L"A(\omega)", title=L"$\eta=%$(standDev),~ L=%$(totalSites+1),~ \varepsilon \in [%$(lowerCutOff),~ %$(bandwidth)],~$%$(disc)",)# yscale=log10)
    ax2 = Axis(f[1,1], width=Relative(0.4), height=Relative(0.4), halign=0.1, valign=0.8, yticklabelsize=24, xticklabelsize=24)
    specFuncIter = IterResults(hamFlow, totalSites, specFuncDefDict, maxSize)
    lines!(ax, freqValues[abs.(freqValues) .< freqValuesZoomLimit], 1e-5 .+ specFuncIter[abs.(freqValues) .< freqValuesZoomLimit], label=L"$J=%$(kondoJ)$", linewidth=4)
    lines!(ax2, freqValues, 1e-5 .+ specFuncIter, linewidth=2)
    axislegend(ax)
    save("specFunc-Kondo-kspace-$(totalSites)-$(disc)-$(lowerCutOff)-$(bandwidth)-$(maxSize).pdf", f)

    return hamFlow, specFuncDefDict
end


function PseudoGapping(
        kondoJ::Float64,
        disc::String,
        totalSites::Int64,
        lowerCutOff::Float64,
        bandwidth::Float64,
        bathInt::Float64,
        cavityRate::Int64,
    )
    @assert (totalSites - initSites) % addPerStep == 0
    specFuncDefDict = Dict("create" => Tuple{String, Vector{Int64}, Float64}[],
                           "destroy" => Tuple{String, Vector{Int64}, Float64}[])
    for site in 1:totalSites
        push!(specFuncDefDict["create"], ("+-+", [2, 1, 2 * site + 1], 0.5))
        push!(specFuncDefDict["create"], ("+-+", [1, 2, 2 * site + 2], 0.5))
        push!(specFuncDefDict["destroy"], ("+--", [1, 2, 2 * site + 1], 0.5))
        push!(specFuncDefDict["destroy"], ("+--", [2, 1, 2 * site + 2], 0.5))
    end
    #=specFuncDefDict["create"] = [("+-", [1,2], 1.)]=#
    #=specFuncDefDict["destroy"] = [("+-", [2, 1], 1.)]=#

    dispersion = zeros(totalSites)
    if disc == "log"
        dispersion[2:2:end] .= 10. .^ range(log10(lowerCutOff), stop=log10(bandwidth), length=div(totalSites-1, 2))
        dispersion[3:2:end] .= -1 .* dispersion[2:2:end]
    else
        dispersion[1:2:end] = range(abs(lowerCutOff), stop=abs(bandwidth), length=div(totalSites+1, 2)) 
        dispersion[2:2:end] .= -1 .* dispersion[3:2:end]
    end
    minceIndices = collect(2 * (1 + initSites):2 * addPerStep:2 * (1 + totalSites))

    f = Figure()
    ax = Axis(f[1,1], xlabel=L"\omega",ylabel=L"A(\omega)", title=L"$\eta=%$(standDev),~ L=%$(totalSites+1),~%$(bathInt),~ \varepsilon \in [%$(lowerCutOff),~ %$(bandwidth)],~$%$(disc)",)# yscale=log10)
    #=ax2 = Axis(f[1,1], width=Relative(0.4), height=Relative(0.4), halign=0.1, valign=0.8, yticklabelsize=24, xticklabelsize=24)=#
    for cavityHead in totalSites+1:-cavityRate:1
        kondoModel = KondoModel(dispersion, kondoJ, bathInt; intLegs=2, globalField=-1e-8, cavityIndices=collect(totalSites:-1:cavityHead))
        hamFlow = MinceHamiltonian(kondoModel, minceIndices)
        specFuncIter = IterResults(hamFlow, totalSites, specFuncDefDict, maxSize)
        scatter!(ax, freqValues[abs.(freqValues) .< freqValuesZoomLimit], 1e-5 .+ specFuncIter[abs.(freqValues) .< freqValuesZoomLimit], label=L"%$(cavityHead)", markersize=8)
    end
    axislegend(ax)
    save("specFunc-Kondo-PG-$(bathInt)-$(totalSites)-$(disc)-$(lowerCutOff)-$(bandwidth)-$(maxSize).pdf", f)
end


function ExactResults(hamFlow, totalSites::Int64, specFuncDefDict)
    totalSpecFunc = zeros(length(freqValues))
    for (i, num) in enumerate(initSites:addPerStep:totalSites)
        if num < totalSites
            continue
        end
        basis = BasisStates(2 * (1 + num); 
                            #=localCriteria=x->x[1]+x[2]==1,=#
                            #=totOccReq=[num, 1 + num, 2 + num]=#
                           )
        fullHam = vcat(hamFlow[1:i]...)
        E, X = Spectrum(fullHam, basis)
        @time specFunc = SpecFunc(E, X,
                            specFuncDefDict,
                            freqValues, basis, standDev, 
                            #=['N'], =#
                            #=(1+num,);=#
                            #=symmetrise=true,=#
                           )
        totalSpecFunc .+= specFunc
    end
    return totalSpecFunc
end

addPerStep = 2
initSites = 1
@assert isodd(initSites)
standDev = 0.1
totalSites = 21
kondoJ = 1e1
bandwidth = 1.
freqValues = collect(-1100 * standDev:10 * standDev:1100 * standDev) .* bandwidth
maxSize = 500
globalField = - 0 * bandwidth/10000

specFuncIter = [[zeros(length(freqValues)), zeros(length(freqValues))] for _ in 1:3]


f = Figure()
ax = Axis(f[1,1], xlabel=L"\omega/D",ylabel=L"A_\mathrm{loc}(\omega)", title=L"$\eta=%$(standDev),~ L=%$(totalSites+1)$")

@assert (totalSites - initSites) % addPerStep == 0
specFuncDefDict = Dict("create" => Tuple{String, Vector{Int64}, Float64}[],
                       "destroy" => Tuple{String, Vector{Int64}, Float64}[])
push!(specFuncDefDict["create"], ("+", [1], 1.))
push!(specFuncDefDict["destroy"], ("-", [1], 1.))
push!(specFuncDefDict["create"], ("+", [2], 1.))
push!(specFuncDefDict["destroy"], ("-", [2], 1.))
#=for site in 1:totalSites=#
#=    push!(specFuncDefDict["create"], ("+-+", [2, 1, 2 * site + 1], 0.5))=#
#=    push!(specFuncDefDict["create"], ("+-+", [1, 2, 2 * site + 2], 0.5))=#
#=    push!(specFuncDefDict["destroy"], ("+--", [1, 2, 2 * site + 1], 0.5))=#
#=    push!(specFuncDefDict["destroy"], ("+--", [2, 1, 2 * site + 2], 0.5))=#
#=end=#

minceIndices = collect(2 * (1 + initSites):2 * addPerStep:2 * (1 + totalSites))

disc = ["log", "lin"]
upperCutOff = [bandwidth, bandwidth/1e3]
lowerCutOff = [1e-6, 0.]
for i in 1:1
    dispersion = Dispersion(totalSites, lowerCutOff[i], upperCutOff[i], disc[i]; phSymmetry=true)
    kondoModel = KondoModel(dispersion, kondoJ; globalField=globalField)
    append!(kondoModel, [("n", [1], -1e0), ("n", [2], -1e0), ("nn", [1, 2], 2e0)])
    hamFlow = MinceHamiltonian(kondoModel, minceIndices)
    specFuncIter[i][1] .= IterResults(hamFlow, totalSites, specFuncDefDict, maxSize)
    kondoModel = KondoModel(dispersion, 0.; globalField=globalField)
    append!(kondoModel, [("n", [1], -1e2), ("n", [2], -1e2), ("nn", [1, 2], 2e2)])
    hamFlow = MinceHamiltonian(kondoModel, minceIndices)
    specFuncIter[i][2] .= IterResults(hamFlow, totalSites, specFuncDefDict, maxSize)
end

labels = ["logarithmic", "quasi-deg."]
linestyles = [:solid, :dash]
for i in 1:2
    scatter!(ax, freqValues ./ bandwidth, 1e-5 .+ specFuncIter[i][1], label=L"%$(labels[i])$$")
    lines!(ax, freqValues ./ bandwidth, 1e-5 .+ specFuncIter[i][2], linestyle=linestyles[i])
end

axislegend(ax, labelsize=24)
save("specFunc-Kondo-kspace-$(totalSites)-$(bandwidth)-$(maxSize).pdf", f)


#=totalSites = 21=#
#=k_points = 1:2:totalSites=#
#=specFuncIter = [zeros(length(freqValues)) for _ in k_points]=#
#==#
#=disc = "log"=#
#=upperCutOff = bandwidth=#
#=lowerCutOff = bandwidth/1e4=#
#==#
#=f = Figure()=#
#=ax = Axis(f[1,1], xlabel=L"\omega/D",ylabel=L"A_\mathrm{d,\vec{k}}(\omega)", title=L"$\eta=%$(standDev),~ L=%$(totalSites+1)$")=#
#==#
#=for (i, point) in enumerate(k_points)=#
#=    specFuncDefDict = Dict("create" => Tuple{String, Vector{Int64}, Float64}[],=#
#=                           "destroy" => Tuple{String, Vector{Int64}, Float64}[])=#
#=    for site in 1:totalSites=#
#=        push!(specFuncDefDict["create"], ("+-+", [2, 1, 2 * site + 1], 0.5)) # S_d^- c^†_up=#
#=        push!(specFuncDefDict["create"], ("+-+", [1, 2, 2 * site + 2], 0.5)) # S_d^+ c^†_down=#
#=        push!(specFuncDefDict["create"], ("+-", [1 + 2 * point, 2 * site + 1], 0.5)) # c^†_kup c_up=#
#=        push!(specFuncDefDict["create"], ("+-", [2 + 2 * point, 2 * site + 2], 0.5)) # c^†_kdown c_down=#
#=    end=#
#=    #=push!(specFuncDefDict["destroy"], ("n-", [2 * totalSites + 1, 1 + 2 * point], 0.5))  # c_{k up}=#=#
#=    #=push!(specFuncDefDict["destroy"], ("h-", [2 * totalSites + 1, 1 + 2 * point], 0.5))  # c_{k up}=#=#
#=    #=push!(specFuncDefDict["destroy"], ("n-", [2 * totalSites + 1, 2 + 2 * point], 0.5))  # c_{k down}=#=#
#=    #=push!(specFuncDefDict["destroy"], ("h-", [2 * totalSites + 1, 2 + 2 * point], 0.5))  # c_{k down}=#=#
#==#
#=    minceIndices = collect(2 * (1 + initSites):2 * addPerStep:2 * (1 + totalSites))=#
#==#
#=    dispersion = Dispersion(totalSites, lowerCutOff, upperCutOff, disc; phSymmetry=true)=#
#=    kondoModel = KondoModel(dispersion, kondoJ; globalField=globalField)=#
#=    hamFlow = MinceHamiltonian(kondoModel, minceIndices)=#
#=    specFuncIter[i] = IterResults(hamFlow, totalSites, specFuncDefDict, maxSize)=#
#=    if sum(specFuncIter[i]) > 1e-2=#
#=        scatter!(ax, freqValues ./ bandwidth, 1e-5 .+ specFuncIter[i], label=L"%$(point)$$")=#
#=    end=#
#=end=#

#=axislegend(ax, labelsize=24)=#
#=save("specFunc-Kondo-kspace-Adk-$(totalSites)-$(bandwidth)-$(maxSize).pdf", f)=#
