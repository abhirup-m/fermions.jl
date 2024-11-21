using fermions, Plots, Measures
include("../src/iterDiag.jl")
include("../src/correlations.jl")

totalSites = 5
initSites = 1
addPerStep = 2
maxSize = 500
bandwidth = 1.
hybridisation = 0.5
U = 5.
Ed = -U/2
freqValues = collect(range(-5. * bandwidth, stop=5. * bandwidth, length=1000))
standDev = 0.1


function IterResults(hamFlow, totalSites::Int64, specFuncDefDict, maxSize::Int64)
    savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, maxSize;
                                                         symmetries=Char['N'],#, 'S'],
                                                         specFuncDefDict=specFuncDefDict,
                                                         #=occReq=(x,N)->abs(x-div(N,2)) ≤ 4=#
                                                         )
    totalSpecFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev;
                                 degenTol=1e-10,
                                 #=occReq=(x,N) -> x == div(N,2),=#
                                 #=excOccReq=(x,N) -> abs(x - div(N,2)) == 1,=#
                                 #=symmetrise=true,=#
                           )
    return totalSpecFunc
end


#=dispersion = zeros(totalSites)=#
#=dispersion[1:2:end] = range(0, stop=bandwidth, length=div(totalSites+1, 2))=#
#=dispersion[2:2:end] = -1 .* dispersion[3:2:end]=#
#=siam = SiamKondoKSpace(dispersion, hybridisation, 0., Ed, U)=#
#=minceIndices = collect(2 * (1 + initSites):2 * addPerStep:2 * (1 + totalSites))=#
#=hamFlow = MinceHamiltonian(siam, minceIndices)=#
#=specFuncDefDict = Dict("create" => Tuple{String, Vector{Int64}, Float64}[],=#
#=                       "destroy" => Tuple{String, Vector{Int64}, Float64}[])=#
#=for site in 3:2:5=#
#=    push!(specFuncDefDict["create"], ("+", [1], 1.))=#
#=    push!(specFuncDefDict["destroy"], ("-", [site], 1.))=#
#==#
#=    specFuncIter = IterResults(hamFlow, totalSites, specFuncDefDict, maxSize)=#
#=    specFuncIter ./= sum(specFuncIter .* (maximum(freqValues) - minimum(freqValues)) ./ (length(freqValues) - 1))=#
#=    println((site, specFuncIter[end]))=#
#=end=#


dispersion = zeros(totalSites)
dispersion[1:2:end] = range(0, stop=bandwidth, length=div(totalSites+1, 2))
dispersion[2:2:end] = -1 .* dispersion[3:2:end]
basisStates = BasisStates(2 * (1 + totalSites))
siam = SiamKondoKSpace(dispersion, hybridisation, 0., Ed, U)
E, X = Spectrum(siam, basisStates)
specFuncDefDict = Dict("create" => Tuple{String, Vector{Int64}, Float64}[],
                       "destroy" => Tuple{String, Vector{Int64}, Float64}[])
push!(specFuncDefDict["create"], ("+", [1], 1.))
push!(specFuncDefDict["destroy"], ("-", [3], 1.))
totalSpecFunc = zeros(length(freqValues))
gstate = X[sortperm(E)[1]]
broadeningFunc(x, standDev) = standDev ./ (x .^ 2 .+ standDev .^ 2)
for (Ei, Xi) in zip(E, X)
    if Ei == minimum(E)
        continue
    end
    overlap1 = StateOverlap(gstate, ApplyOperator(specFuncDefDict["destroy"], Xi)) * StateOverlap(Xi, ApplyOperator(specFuncDefDict["create"], gstate))
    overlap2 = StateOverlap(gstate, ApplyOperator(specFuncDefDict["create"], Xi)) * StateOverlap(Xi, ApplyOperator(specFuncDefDict["destroy"], gstate))
    if abs(overlap2) > 1e-5 #|| abs(overlap1) > 1e-10
        println(overlap2)
    end
end
