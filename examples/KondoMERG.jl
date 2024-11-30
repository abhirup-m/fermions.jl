using fermions

function RGFlow(
        energies::Vector{Float64},
        kondoJ::Float64,
        dofStates::Function,
    )
    @assert length(energies) > 1
    @assert issorted(energies, rev=true)
    kondoJFlow = [kondoJ]
    for (i, energy) in enumerate(energies)
        if energy - kondoJFlow[end]/4 < 0
            return kondoJFlow, i-1
        end
        if kondoJFlow[end] < 0
            kondoJFlow[end] = 0
            return kondoJFlow, i-1
        end
        deltaEnergy = i < length(energies) ? energies[i] - energies[i+1] : energies[i-1] - energies[i]
        deltaKondoJ = dofStates(energy, maximum(energies)) * kondoJFlow[end]^2 * deltaEnergy / (energy - kondoJFlow[end]/4)
        push!(kondoJFlow, deltaKondoJ + kondoJFlow[end])
    end
    return kondoJFlow, length(energies)
end

D = 3.
deltaD = 0.01
kondoJ = 0.2 * D
dofStates = (E,D) -> (1/D)*(1 - E^2/D^2)^0.5
energies = collect(D:-deltaD:0)
kondoJFlow, fixedPointStep = RGFlow(energies, kondoJ, dofStates)

totalSites = 11
initStates = 5
kondoJMatrix = zeros(totalSites, totalSites)
kondoJMatrix[1:initStates,1:initStates] .= kondoJFlow[end]
for i in 1:(totalSites - initStates)
    kondoJMatrix[initStates+i, 1:initStates+i] .= kondoJFlow[end - i]
    kondoJMatrix[1:initStates+i, initStates+i] .= kondoJFlow[end - i]
end
hamltonian = KondoModel(energies[end:-1:end-totalSites+1], kondoJMatrix)
hamiltonianFlow = MinceHamiltonian(hamltonian, collect(2 + 2 * initStates:4:2 + 2 * totalSites))
savePaths, resultsDict = IterDiag(hamiltonianFlow, 500;
                     symmetries=Char['N', 'S'],
                     correlationDefDict=Dict("SF-d0" => [("+-+-", [1,2,4,3], 1.), ("+-+-", [2,1,3,4], 1.)]),
                     calculateThroughout=true,
                    ) 
display(resultsDict)
