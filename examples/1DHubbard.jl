#### Iterative diagonalisation solution of the 1D Hubbard model ####
using ProgressMeter, LinearAlgebra, Plots, Measures
include("../src/base.jl")
include("../src/correlations.jl")
include("../src/iterativeDiag.jl")

function dimerHamiltonian(U, t)
    hopping = Dict(
                   ("+-", [1, 3]) => -t,
                   ("+-", [2, 4]) => -t,
                   ("+-", [3, 1]) => -t,
                   ("+-", [4, 2]) => -t,
                  )
    correlation = Dict(
                       ("n", [1]) => -U/2,
                       ("n", [2]) => -U/2,
                       ("nn", [1, 2]) => U,
                       ("n", [3]) => -U/2,
                       ("n", [4]) => -U/2,
                       ("nn", [3, 4]) => U,
                  )
    return mergewith(+, hopping, correlation)
end


function dimerAdditionalHamiltonian(U, t, index)
    upPosition = 2 * index - 1
    hopping = Dict(
                   ("+-", [upPosition, upPosition - 2]) => -t,
                   ("+-", [upPosition + 1, upPosition - 1]) => -t,
                   ("+-", [upPosition - 2, upPosition]) => -t,
                   ("+-", [upPosition - 1, upPosition + 1]) => -t,
                  )
    
    correlation = Dict(
                       ("n", [upPosition]) => -U/2,
                       ("n", [upPosition + 1]) => -U/2,
                       ("nn", [upPosition, upPosition + 1]) => U
                      )
    return mergewith(+, hopping, correlation)
end

function main(numSteps, U)
    t = 1.0
    initBasis = BasisStates(4; occCriteria=x -> x ∈ [1, 2, 3], magzCriteria = x -> x ∈ [-1, 0, 1])#; magzCriteria=magzCriteria, occCriteria=occCriteria)
    retainSize = 50
    hamiltonianFamily = [dimerHamiltonian(U, t)]
    numStatesFamily = Int64[2]
    for i in 1:numSteps
        push!(hamiltonianFamily, dimerAdditionalHamiltonian(U, t, 2 + i)
             )
        push!(numStatesFamily, 2 + i)
    end
    spectrumFamilyGS = iterativeDiagonaliser(hamiltonianFamily, initBasis, numStatesFamily, retainSize;
                                           #occCriteria= (x,y) -> x ∈ [y-1, y, y+1], magzCriteria = x -> x ∈ [-1, 0, 1]
                                          )
    return
    # spectrumFamilyES = iterativeDiagonaliser(hamiltonianFamily, initBasis, numStatesFamily, retainSize;
                                           # occCriteria= (x,y) -> x ∈ [y-1], magzCriteria = x -> x ∈ [-2, -1, 0], keepfrom="UV"
                                          # )
    # spectrumFamilyESDag = iterativeDiagonaliser(hamiltonianFamily, initBasis, numStatesFamily, retainSize;
                                           # occCriteria= (x,y) -> x ∈ [y+1], magzCriteria = x -> x ∈ [0, 1, 2], keepfrom="UV"
                                          # )
    plots = []
    freqArray = collect(range(-10.0, stop=10.0, length=500))
    specfunc = 0 .* freqArray
    broadening = 0.01
    names = []
    for (numStates, (energyGS, statesGS)) in zip(numStatesFamily, spectrumFamilyGS)
        if numStates % 2 == 1
            continue
        end
        sector = (numStates, 0)
        probe = Dict(("-", [1]) => 1.0)
        probeDag = Dict(("+", [1]) => 1.0)
        specfunc += specFunc(statesGS[sector][1], energyGS[sector][1], (energyGS, energyGS), (statesGS, statesGS), probe, probeDag, freqArray, broadening)
        p = plot(freqArray, specfunc, linewidth=2, title="\$N=$numStates,\\quad U/t=$(U/t)\$", xlabel="\$\\omega/t\$", ylabel="\$A_1(\\omega)\$", size=(450, 200), leftmargin=2mm, legend=false)
        savefig(p, "specfunc_$(numStates).pdf")
        push!(names, "specfunc_$(numStates).pdf")
    end
    run(`pdfunite $names specfunc.pdf`)
    run(`rm $names`)
end

main(5, 0.0)
