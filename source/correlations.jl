using LinearAlgebra
include("fermionise.jl")

function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    minimumEnergies = minimum.(values(eigvals))
    minimumBlock = collect(keys(eigvals))[argmin(minimumEnergies)]
    println(collect(keys(eigvals)), minimumBlock)
    println(minimumEnergies)
    minimumIndex = argmin(eigvals[minimumBlock])
    gstate = eigvecs[minimumBlock][minimumIndex]

    for (state, coeff) in zip(basisStates[minimumBlock], gstate)
        if round(coeff, digits=8) == 0
            continue
        end
        println(state, round(coeff, digits=8))
    end

    spinFlipCorrOplist = Dict(("n", [4]) => 0.5, ("n", [4]) => 0.5)
    operatorMatrix = generalOperatorMatrix(basisStates, spinFlipCorrOplist)
    correlation = gstate' * operatorMatrix[minimumBlock] * gstate
    spinFlipCorrOplist = Dict(("n", [6]) => 0.5, ("n", [6]) => 0.5)
    operatorMatrix = generalOperatorMatrix(basisStates, spinFlipCorrOplist)
    newcorrelation = gstate' * operatorMatrix[minimumBlock] * gstate
    spinFlipCorrOplist = Dict(("n", [8]) => 0.5, ("n", [8]) => 0.5)
    operatorMatrix = generalOperatorMatrix(basisStates, spinFlipCorrOplist)
    newcorrelation2 = gstate' * operatorMatrix[minimumBlock] * gstate
    println("----------------- ($correlation, $newcorrelation, $newcorrelation2) -----------------")
    return correlation, gstate
end

function gstateCorrelation(basisStates::Dict{Tuple{Int64,Int64},Vector{BitArray}}, eigvals::Dict{Tuple{Int64,Int64},Vector{Float64}}, eigvecs::Dict{Tuple{Int64,Int64},Vector{Vector{Float64}}}, correlationOperatorArr::Vector{Dict{Tuple{String,Vector{Int64}},Float64}})
    minimumEnergies = minimum.(values(eigvals))
    println("ME", minimumEnergies, collect(keys(eigvals)))
    minimumBlock = collect(keys(eigvals))[argmin(minimumEnergies)]
    correlationArr = []
    for minimumBlock in [(4, 2), (4, -2), (4, 0)]
        minimumIndex = argmin(eigvals[minimumBlock])
        println(minimumBlock, eigvals[minimumBlock])
        gstate = eigvecs[minimumBlock][minimumIndex]
        gstate2 = eigvecs[minimumBlock][minimumIndex+1]
        for (state, coeff) in zip(basisStates[minimumBlock], gstate2)
            if abs(coeff) > 1e-5
                println((prettyPrint(state), coeff))
            end
        end
        println("--------")
        correlationArr = [gstate2' * generalOperatorMatrix(basisStates, correlationOperator)[minimumBlock] * gstate2
                          for correlationOperator in correlationOperatorArr]
        println("2", correlationArr)
        correlationArr = [gstate' * generalOperatorMatrix(basisStates, correlationOperator)[minimumBlock] * gstate
                          for correlationOperator in correlationOperatorArr]
        println("1", correlationArr)
        println("--------")
    end
    return round.(correlationArr, digits=trunc(Int, -log10(TOLERANCE)))
end


function simpleCorrelation(basis, eigvals, eigvecs, corrDef)
    gsEnergy = minimum(round.(minimum.(values(eigvals)), digits=trunc(Int, -log10(TOLERANCE))))
    corrOp = generalOperatorMatrix(basis, corrDef)
    correlations = []
    for (block, eigvalsBlock) in eigvals
        for (i, _) in enumerate(eigvalsBlock[abs.(eigvalsBlock .- gsEnergy).<TOLERANCE])
            correlation = eigvecs[block][i]' * corrOp[block] * eigvecs[block][i]
            push!(correlations, correlation)
        end
    end
    return correlations
end
