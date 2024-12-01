using fermions, CairoMakie, Measures
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

mergResults = Dict("SF-d0" => Float64[])
initStates = 1
for totalSites in [3, 21]
    kondoJMatrix = zeros(1 + 2 * totalSites, 1 + 2 * totalSites)
    kondoJMatrix[1:1 + 2 * initStates, 1:1 + 2 * initStates] .= kondoJFlow[end]
    for i in 1:totalSites - initStates
        kondoJMatrix[-1 + 2 * initStates + 2 * i:1 + 2 * initStates + 2 * i, 1:1 + 2 * initStates + 2 * i] .= kondoJFlow[end - i]
        kondoJMatrix[1:1 + 2 * initStates + 2 * i, -1 + 2 * initStates + 2 * i:1 + 2 * initStates + 2 * i] .= kondoJFlow[end - i]
    end
    println(kondoJMatrix)
    println([[energies[end]]; repeat(energies[end-1:-1:end-totalSites], inner=2)])
    hamltonian = KondoModel([[energies[end]]; repeat(energies[end-1:end-totalSites], inner=2)], kondoJMatrix)
    hamiltonianFlow = MinceHamiltonian(hamltonian, collect((4 + 2 * initStates):2:(4 + 2 * totalSites)))
    savePaths, resultsDict = IterDiag(hamiltonianFlow, 500;
                         symmetries=Char['N', 'S'],
                         correlationDefDict=Dict("SF-d0" => [("+-+-", [1,2,4,3], 1.), ("+-+-", [2,1,3,4], 1.)]),
                         #=calculateThroughout=true,=#
                        ) 
    push!(mergResults["SF-d0"], resultsDict["SF-d0"])
end
p = CairoMakie.plot(mergResults["SF-d0"])
save("merg.pdf", p)
