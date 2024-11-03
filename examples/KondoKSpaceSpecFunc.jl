using fermions, CairoMakie, Measures, ProgressMeter
include("../src/iterDiag.jl")

set_theme!(merge(theme_ggplot2(), theme_latexfonts()))
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


function IterResults(hamFlow, totalSites::Int64, specFuncDefDict, maxSize::Int64,)
    savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                                                         symmetries=Char['N'],#, 'S'],
                                 specFuncDefDict=specFuncDefDict,
                                 #=occReq=(x,N)->abs(x-div(N,2)) ≤ 4=#
                                )
    totalSpecFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                                       #=occReq=(x,N) -> x == div(N,2),=#
                                       #=excOccReq=(x,N) -> abs(x - div(N,2)) == 1,=#
                                       #=symmetrise=true,=#
                           )
    return totalSpecFunc
end


function LargerSystem(
        kondoJVals::Vector{Float64},
        disc::String,
        totalSitesVals::Vector{Int64},
        lowerCutOff::Float64,
        bandwidth::Float64,
    )
    for totalSites in totalSitesVals
        @assert (totalSites - initSites) % addPerStep == 0
        specFuncDefDict = Dict("create" => Tuple{String, Vector{Int64}, Float64}[],
                               "destroy" => Tuple{String, Vector{Int64}, Float64}[])
        for site in 1:totalSites
            push!(specFuncDefDict["create"], ("+-+", [2, 1, 2 * site + 1], 0.5))
            push!(specFuncDefDict["create"], ("+-+", [1, 2, 2 * site + 2], 0.5))
            push!(specFuncDefDict["destroy"], ("+--", [1, 2, 2 * site + 1], 0.5))
            push!(specFuncDefDict["destroy"], ("+--", [2, 1, 2 * site + 2], 0.5))
        end

        dispersion = zeros(totalSites)
        if disc == "log"
            dispersion[2:2:end] .= 10. .^ range(log10(lowerCutOff), stop=log10(bandwidth), length=div(totalSites-1, 2))
            dispersion[3:2:end] .= -1 .* dispersion[2:2:end]
        else
            dispersion[1:2:end] = range(abs(lowerCutOff), stop=abs(bandwidth), length=div(totalSites+1, 2)) 
            dispersion[2:2:end] .= -1 .* dispersion[3:2:end]
        end
        println((dispersion[1], dispersion[end]))
        minceIndices = collect(2 * (1 + initSites):2 * addPerStep:2 * (1 + totalSites))

        f = Figure()
        ax = Axis(f[1,1], xlabel=L"\omega",ylabel=L"A(\omega)", title=L"$\eta=%$(standDev),~ L=%$(totalSites+1),~ \varepsilon \in [%$(lowerCutOff),~ %$(bandwidth)],~$%$(disc)",)# yscale=log10)
        ax2 = Axis(f[1,1], width=Relative(0.4), height=Relative(0.4), halign=0.1, valign=0.8, yticklabelsize=24, xticklabelsize=24)
        for kondoJ in kondoJVals
            kondoModel = KondoModel(dispersion, kondoJ, globalField=-0.)
            hamFlow = MinceHamiltonian(kondoModel, minceIndices)
            specFuncIter = IterResults(hamFlow, totalSites, specFuncDefDict, maxSize)
            scatter!(ax, freqValues[abs.(freqValues) .< freqValuesZoomLimit], 1e-5 .+ specFuncIter[abs.(freqValues) .< freqValuesZoomLimit], label=L"$J=%$(kondoJ)$", markersize=8)
            lines!(ax2, freqValues, 1e-5 .+ specFuncIter, linewidth=2)
        end
        axislegend(ax)
        save("specFunc-Kondo-kspace-$(totalSites)-$(disc)-$(lowerCutOff)-$(bandwidth)-$(maxSize).pdf", f)
    end
end


addPerStep = 2
initSites = 1
@assert isodd(initSites)
standDev = 0.02
freqValues = collect(-2:0.005:2)
freqValuesZoomLimit = 0.3

maxSize = 1000
@sync begin
@async LargerSystem([1e2, 1e0, 1e-2, 1e-4, 0.], "lin", [21], 0., 0.0)
@async LargerSystem([1e2, 1e0, 1e-2, 1e-4, 0.], "lin", [21], 0., 0.1)
@async LargerSystem([1e2, 1e0, 1e-2, 1e-4, 0.], "log", [21], 1e-4, 0.1)
@async LargerSystem([1e2, 1e0, 1e-2, 1e-4, 0.], "log", [21], 1e-6, 0.1)
end
