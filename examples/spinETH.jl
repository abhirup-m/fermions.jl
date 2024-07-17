using Plots, Measures, ProgressMeter
include("../src/base.jl")
include("../src/correlations.jl")
include("../src/eigen.jl")

theme(:dark)
default(linewidth=3, thickness_scaling=1.5)

timeMax = 20
deltaTime = 0.2
timeValues = range(0, stop=timeMax, step=deltaTime)
probe = [("n", [1], 0.5), ("h", [1], -0.5)]
numSitesSet = [4, 8, 12]
saveNames = []
for numSites in numSitesSet
    @assert numSites % 2 == 0
    basisStates = BasisStates(numSites; magzCriteria=x->sum(x) == numSites/2)
    probeMatEle = [ifelse(collect(keys(state))[1] == 1, 0.5, -0.5) for state in basisStates]

    initState1 = Dict(BitVector(vcat(repeat([[1, 0]], outer=div(numSites, 2))...)) => 1.0)
    initState2 = Dict(BitVector(vcat(repeat([[0, 1]], outer=div(numSites, 2))...)) => 1.0)

    state1Coeffs = fetch.([Threads.@spawn StateOverlap(bstate, initState1) for bstate in basisStates])
    state2Coeffs = fetch.([Threads.@spawn StateOverlap(bstate, initState2) for bstate in basisStates])

    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for i in 1:(numSites-1)
        push!(hamiltonian,("+-", [i + 1, i], 0.5))
        push!(hamiltonian,("+-", [i, i + 1], 0.5))
    end
    push!(hamiltonian,("+-", [numSites, 1], 0.5))
    push!(hamiltonian,("+-", [1, numSites], 0.5))

    hamMatrix = OperatorMatrix(basisStates, hamiltonian)
    probeMatrix = OperatorMatrix(basisStates, probe)

    getOperator(i) = ((I - (hamMatrix .* 0.5im .* deltaTime)) / (I + (hamMatrix .* 0.5im .* deltaTime)))^i
    unitaryOperators = @showprogress desc="N=$(numSites), unitaries" [getOperator(i) for (i, _) in enumerate(timeValues)]

    probeResults = [zeros(length(timeValues)), zeros(length(timeValues))]
    matrixMul1(operator) = state1Coeffs' * (operator * probeMatrix * operator') * state1Coeffs
    matrixMul2(operator) = state2Coeffs' * (operator * probeMatrix * operator') * state2Coeffs
    probeResults[1] = @showprogress desc="N=$(numSites), first state" [matrixMul1(operator) for operator in unitaryOperators]
    probeResults[2] = @showprogress desc="N=$(numSites), second state" [matrixMul2(operator) for operator in unitaryOperators]

    p = plot(range(0, stop=timeMax, step=deltaTime), probeResults; legend=:topright, title="\$N=$(numSites)\$",
             labels=["\$\\uparrow\\uparrow\\ldots\\downarrow\\downarrow\$" "\$\\downarrow\\downarrow\\ldots\\uparrow\\uparrow\$"])
    savefig(p, "spinThermal_$(numSites).pdf")
    push!(saveNames, "spinThermal_$(numSites).pdf")
end
run(`pdfunite $saveNames spinThermal.pdf`)
run(`rm $saveNames`)
