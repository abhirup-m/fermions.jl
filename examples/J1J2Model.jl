using fermions, Serialization
include("../src/modelHamiltonians.jl")

totalSites = 33
@assert totalSites % 3 == 0
partitions = collect(3:3:totalSites)
J1byJ2Values = 0:0.1:2
gapValues = zeros(length(J1byJ2Values))
@time for (i, J1byJ2) in enumerate(J1byJ2Values)
    hamiltonian = J1J2Model(J1byJ2, totalSites)
    hamiltonianFlow = MinceHamiltonian(hamiltonian, partitions)
    savePaths, results = IterDiag(hamiltonianFlow, 2000)
    eigVals = deserialize(savePaths[end-1])["eigVals"]
    gapValues[i] = eigVals[2] - eigVals[1]
end
display(gapValues)
