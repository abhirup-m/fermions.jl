function KondoHamiltonian(basisStates, numBathSites, couplings)
    if length(basisStates) == 0
        println("Chosen basis is empty, bailing out.")
        return
    end
	(kineticEnergy, kondoCoupling) = couplings
	@assert length(kineticEnergy) == numBathSites
	kineticEnergy = repeat(kineticEnergy, inner=2)
	upIndices = range(3, 2 * numBathSites + 1, step=2)

	operatorList = []
	
	for (k1, k2) in Iterators.product(upIndices, upIndices)
		push!(operatorList, (0.25 * kondoCoupling / numBathSites, "n+-", (1, k1, k2)))
		push!(operatorList, (-0.25 * kondoCoupling / numBathSites, "n+-", (1, k1 + 1, k2 + 1)))
		push!(operatorList, (-0.25 * kondoCoupling / numBathSites, "n+-", (2, k1, k2)))
		push!(operatorList, (0.25 * kondoCoupling / numBathSites, "n+-", (2, k1 + 1, k2 + 1)))
		push!(operatorList, (0.5 * kondoCoupling / numBathSites, "+-+-", (1, 2, k1 + 1, k2)))
		push!(operatorList, (0.5 * kondoCoupling / numBathSites, "+-+-", (2, 1, k1, k2 + 1)))
	end

    for upIndex in upIndices
		push!(operatorList, (kineticEnergy[upIndex - 2], "n", (upIndex,)))
		push!(operatorList, (kineticEnergy[upIndex - 1], "n", (upIndex + 1,)))
	end

	return GeneralOperatorMatrix(basisStates, operatorList)
end


function TOSHamiltonianRSpace(basisStates, couplings, show_progress=true)
    # imup    imdn    A1up    A1dn    B1up    B1dn    A2up    A2dn    B2up    B2dn
    # 0       1       2       3       4       5       6       7       8       9
    # forward hoppings: 2 --> 6, 3 --> 7, 4 --> 8, 5 --> 9 
    # A-->B hoppings: 2 --> 4, 3 --> 5, 6 --> 8, 7 --> 9
    # hybridisation: 0 --> 2, 0 --> 4, 1 --> 3, 1 --> 5
    # impurity correlation: n: 0, n: 1, nn: 01
    # bath correlation: n: 2, n: 3, nn: 23, n: 4, n: 5, nn: 45
    
    t_hop, hop_strength, imp_U, t_perp, bath_Ub, Bfield = couplings
    numSitesChain = div(length(basisStates[1]) - 2, 4)
    operatorList = []

    if Bfield != 0 
    for i in 1:2 + 4 * numSitesChain 
        push!(operatorList, [(-1)^i * Bfield, "n", (i,)])
    end
    end

    
    if t_hop != 0 
    for i in 3:6 + 4 * (numSitesChain - 2)
        push!(operatorList, [-t_hop, "+-", (i, i+4)])
        push!(operatorList, [-t_hop, "+-", (i+4, i)])
    end
    end

    if t_perp != 0
    for i in 3:4:3 + 4 * (numSitesChain - 1)
        push!(operatorList, [-t_perp, "+-", (i, i+2)])
        push!(operatorList, [-t_perp, "+-", (i+2, i)])
        push!(operatorList, [-t_perp, "+-", (i+1, i+3)])
        push!(operatorList, [-t_perp, "+-", (i+3, i+1)])
    end
    end

    if hop_strength != 0
    for i in [1 2]
        for j in [2, 4]
            push!(operatorList, [hop_strength, "+-", (i, i + j)])
            push!(operatorList, [hop_strength, "+-", (i + j, i)])
        end
    end
    end

    if imp_U  != 0 || bath_Ub != 0
    for i in [1 3 5]
        corr = i == 1 ? imp_U : bath_Ub
        push!(operatorList, [-0.5 * corr, "n", (i,)])
        push!(operatorList, [-0.5 * corr, "n", (i+1,)])
        push!(operatorList, [corr, "nn", (i, i+1)])
    end
    end
    
    return GeneralOperatorMatrix(basisStates, operatorList, show_progress)
end