using fermions, CairoMakie, Measures, ProgressMeter
#=include("../src/base.jl")=#
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


function IterResults(hamFlow, totalSites::Int64)
    savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                                                         symmetries=Char['N'],#, 'S'],
                                 specFuncDefDict=specFuncDefDict,
                                 occReq=(x,N)->abs(x-div(N,2)) ≤ 4
                                )
    totalSpecFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                                       #=occReq=(x,N) -> x == div(N,2),=#
                                       #=excOccReq=(x,N) -> abs(x - div(N,2)) == 1,=#
                                       #=symmetrise=true,=#
                           )
    return totalSpecFunc
end

function ExactResults(hamFlow, totalSites::Int64)
    totalSpecFunc = zeros(length(freqValues))
    for (i, num) in enumerate(initSites:addPerStep:totalSites)
        basis = BasisStates(2 * (1 + num); 
                            localCriteria=x->x[1]+x[2]==1,
                            totOccReq=[num, 1 + num, 2 + num]
                           )
        fullHam = vcat(hamFlow[1:i]...)
        E, X = Spectrum(fullHam, basis)
        specFunc = SpecFunc(E, X,
                            specFuncDefDict,
                            freqValues, basis, standDev, 
                            ['N'], 
                            (1+num,);
                            #=symmetrise=true,=#
                           )
        totalSpecFunc .+= specFunc
    end
    return totalSpecFunc
end

function BenchMark(
        kondoJVals::Vector{Float64},
        hop_t::Float64,
    )
    f = Figure(size=(700, 800))
    axes = [Axis(f[i, 1], xlabel=L"\omega",ylabel=L"A(\omega)") for i in 1:2]
    errorAxes = [Axis(f[i, 1], width=Relative(0.4), height=Relative(0.4), halign=0.2, valign=0.85, yticklabelsize=20, xticklabelsvisible=false) for i in 1:2]
    axisIndex = 1
    for totalSites in [5, 6]
        for kondoJ in kondoJVals
            kondoModel = KondoModel(totalSites, hop_t, kondoJ, globalField=1e-5)
            hamFlow = MinceHamiltonian(kondoModel, collect(2 * (1 + initSites):2 * addPerStep:2 * (1 + totalSites)))
            specFuncIter = IterResults(hamFlow, totalSites)
            specFuncExact = ExactResults(hamFlow, totalSites)
            lines!(axes[axisIndex], freqValues[freqValues .≥ 0], specFuncIter[freqValues .≥ 0], label=L"ID $J=%$(kondoJ)$")
            scatter!(axes[axisIndex], freqValues[freqValues .≥ 0], specFuncExact[freqValues .≥ 0], linestyle=:dash, label=L"ED $J=%$(kondoJ)$")
            relError = abs.(specFuncIter .- specFuncExact) ./ specFuncExact
            lines!(errorAxes[axisIndex], relError)
        end
        axislegend(axes[axisIndex])
        axisIndex += 1
    end
    save("specFuncComparison.pdf", f)
end


function LargerSystem(
        kondoJVals::Vector{Float64},
        hop_t::Float64
    )
    f = Figure()
    totalSites = 29
    ax = Axis(f[1,1], xlabel=L"\omega",ylabel=L"A(\omega)", title=L"\eta=%$(standDev), L=%$(totalSites+1), t=%$(hop_t)", yscale=log10)
    for kondoJ in kondoJVals
        kondoModel = KondoModel(totalSites, hop_t, kondoJ, globalField=-1e-5)
        hamFlow = MinceHamiltonian(kondoModel, collect(2 * (1 + initSites):2 * addPerStep:2 * (1 + totalSites)))
        specFuncIter = IterResults(hamFlow, totalSites)
        specFuncIter ./= sum(specFuncIter .* (maximum(freqValues) - minimum(freqValues[1])) / (length(freqValues) - 1))
        scatterlines!(ax, freqValues, 1e-5 .+ specFuncIter, label=L"$J=%$(kondoJ)$")
    end
    axislegend(ax)
    save("specFunc-Kondo-real.pdf", f)
end

#=specFuncDefDict = Dict("create" => [("+-+", [1,2,4], 0.5)], "destroy" => [("+--", [2,1,4], 0.5)])=#
specFuncDefDict = Dict("create" => [("+-+", [2,1,3], 0.5), ("+-+", [1,2,4], 0.5)], "destroy" => [("+--", [1,2,3], 0.5), ("+--", [2,1,4], 0.5)])
initSites = 1
maxSize = 1500
hop_t = 0.01
addPerStep = 1
standDev = 0.02
freqValues = collect(0:0.01:2.)
BenchMark([1., 0.], hop_t)

freqValues = collect(-2:0.005:2.)
LargerSystem([1., 0.5, 0.], hop_t)
