#### Iterative diagonalisation solution of the 1D Hubbard model ####
using ProgressMeter, LinearAlgebra, Plots
using fermions

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
    initBasis = fermions.BasisStates(4; occCriteria=x -> x ∈ [1, 2, 3], magzCriteria = x -> x ∈ [-1, 0, 1])#; magzCriteria=magzCriteria, occCriteria=occCriteria)
    retainSize = 50
    hamiltonianFamily = [dimerHamiltonian(U, t)]
    numStatesFamily = Int64[2]
    for i in 1:numSteps
        push!(hamiltonianFamily, dimerAdditionalHamiltonian(U, t, 2 + i)
             )
        push!(numStatesFamily, 2 + i)
    end
    spectrumFamily = fermions.iterativeDiagonaliser(hamiltonianFamily, initBasis, numStatesFamily, retainSize;
                                           occCriteria= (x,y) -> x ∈ [y-1, y, y+1], magzCriteria = x -> x ∈ [-1, 0, 1]
                                          )#; occCriteria=occCriteria, magzCriteria = magzCriteria)
    plots = []
    for (numStates, (E, X)) in zip(numStatesFamily, spectrumFamily)
        if numStates % 2 == 1
            continue
        end
        sector = (numStates, 0)
        groundState = X[sector][1]
        energyGs = E[sector][1]
        probe = Dict(("-", [2 * numStates - 1]) => 1.0)
        probeDag = Dict(("+", [2 * numStates - 1]) => 1.0)
        freqArray = collect(range(-10, stop=10, length=500))
        specfunc = fermions.specFunc(groundState, energyGs, E, X, probe, probeDag, freqArray, 0.0001)

        push!(plots, plot(freqArray, specfunc, linewidth=2, xlabel="\$\\omega\$", ylabel="\$A(\\omega)\$", size=(300, 200 * length(plots))))
    end
    p = plot(plots..., layout=(length(plots), 1))
    savefig(p, "specfunc.pdf")
end

main(10, 1.0)
