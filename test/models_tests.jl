@testset "Kondo model, real-space" begin
    globalField = 1e-5
    @testset for hop_t in [0., -1., 1e-14, 1e-17]
        @testset for kondoJ in [0., 1., 1e-14, 1e-17]
            @testset for bathInt in [1., -1., 1e-17]
                @testset for numBathSites in [1, 2, 10, 19]
                    hamiltonian = KondoModel(numBathSites, hop_t, kondoJ, bathInt; globalField=globalField)
                    comparisonHamiltonian = []
                    for i in 0:numBathSites
                        push!(comparisonHamiltonian, ("n", [1 + 2 * i], globalField))
                        push!(comparisonHamiltonian, ("n", [2 + 2 * i], -globalField))
                    end
                    if abs(hop_t) > 1e-15
                        for i in 1:numBathSites-1
                            push!(comparisonHamiltonian, ("+-", [1 + 2 * i, 3 + 2 * i], -hop_t))
                            push!(comparisonHamiltonian, ("+-", [2 + 2 * i, 4 + 2 * i], -hop_t))
                            push!(comparisonHamiltonian, ("+-", [3 + 2 * i, 1 + 2 * i], -hop_t))
                            push!(comparisonHamiltonian, ("+-", [4 + 2 * i, 2 + 2 * i], -hop_t))
                        end
                    end
                    if abs(kondoJ) > 1e-15
                        push!(comparisonHamiltonian, ("+-+-", [1, 2, 4, 3], kondoJ/2))
                        push!(comparisonHamiltonian, ("+-+-", [2, 1, 3, 4], kondoJ/2))
                        push!(comparisonHamiltonian, ("nn", [1, 3], kondoJ/4))
                        push!(comparisonHamiltonian, ("nn", [1, 4], -kondoJ/4))
                        push!(comparisonHamiltonian, ("nn", [2, 3], -kondoJ/4))
                        push!(comparisonHamiltonian, ("nn", [2, 4], kondoJ/4))
                    end
                    if abs(bathInt) > 1e-15
                        push!(comparisonHamiltonian, ("n", [3], -bathInt/2))
                        push!(comparisonHamiltonian, ("n", [4], -bathInt/2))
                        push!(comparisonHamiltonian, ("nn", [3, 4], bathInt))
                    end
                    @test issetequal(hamiltonian, comparisonHamiltonian)
                end
            end
        end
    end
end

@testset "Kondo model, k-space" begin
    globalField = 1e-5
    @testset for kondoJ in [0., 1., 1e-14, 1e-17]
        @testset for bathInt in [1., -1., 1e-17]
            @testset for numBathSites in [1, 2, 5]
                @testset for bathIntLegs in 1:4
                    dispersion = rand(numBathSites)
                    if numBathSites ≥ 2
                        cavityIndices = [numBathSites-1, numBathSites-3]
                    else
                        cavityIndices = Int64[]
                    end
                    hamiltonian = KondoModel(dispersion, kondoJ, bathInt; bathIntLegs=bathIntLegs, globalField=globalField, cavityIndices=cavityIndices)
                    comparisonHamiltonian = []
                    for i in 0:numBathSites
                        push!(comparisonHamiltonian, ("n", [1 + 2 * i], globalField))
                        push!(comparisonHamiltonian, ("n", [2 + 2 * i], -globalField))
                    end
                    for i in 1:numBathSites
                        if abs(dispersion[i]) > 1e-15
                            push!(comparisonHamiltonian, ("n", [1 + 2 * i], dispersion[i]))
                            push!(comparisonHamiltonian, ("n", [2 + 2 * i], dispersion[i]))
                        end
                    end
                    if abs(kondoJ) > 1e-15
                        for (i,j) in Iterators.product(1:numBathSites, 1:numBathSites)
                            if i ∈ cavityIndices || j ∈ cavityIndices
                                continue
                            end
                            push!(comparisonHamiltonian, ("+-+-", [1, 2, 2 + 2 * i, 1 + 2 * j], kondoJ/2))
                            push!(comparisonHamiltonian, ("+-+-", [2, 1, 1 + 2 * i, 2 + 2 * j], kondoJ/2))
                            push!(comparisonHamiltonian, ("n+-", [1, 1 + 2 * i, 1 + 2 * j], kondoJ/4))
                            push!(comparisonHamiltonian, ("n+-", [1, 2 + 2 * i, 2 + 2 * j], -kondoJ/4))
                            push!(comparisonHamiltonian, ("n+-", [2, 1 + 2 * i, 1 + 2 * j], -kondoJ/4))
                            push!(comparisonHamiltonian, ("n+-", [2, 2 + 2 * i, 2 + 2 * j], kondoJ/4))
                        end
                    end
                    if abs(bathInt) > 1e-15
                        for (i,j, k, l) in Iterators.product(repeat([1:numBathSites], 4)...)
                            if length(unique([i, j, k, l])) > bathIntLegs
                                continue
                            end
                            push!(comparisonHamiltonian, ("+-+-", [2 * i + 1, 2 * j + 1, 2 * k + 1, 2 * l + 1], -bathInt/2))
                            push!(comparisonHamiltonian, ("+-+-", [2 * i + 2, 2 * j + 2, 2 * k + 2, 2 * l + 2], -bathInt/2))
                            push!(comparisonHamiltonian, ("+-+-", [2 * i + 1, 2 * j + 1, 2 * k + 2, 2 * l + 2], bathInt/2))
                            push!(comparisonHamiltonian, ("+-+-", [2 * i + 2, 2 * j + 2, 2 * k + 1, 2 * l + 1], bathInt/2))
                        end
                    end

                    @test issetequal(hamiltonian, comparisonHamiltonian)
                end
            end
        end
    end
end

@testset "Kondo model, k-space, anisotropic" begin
    globalField = 1e-5
    @testset for kondoJ in [0., 1., 1e-14, 1e-17]
        @testset for numBathSites in [1, 2, 5]
            @testset for bathInt in [1., -1., 1e-17]
                @testset for bathIntLegs in 1:4
                    dispersion = rand(numBathSites)
                    kondoJMatrix = rand(numBathSites, numBathSites) .* kondoJ
                    bathIntMatrix = rand(numBathSites, numBathSites, numBathSites, numBathSites) .* bathInt
                    bathIntFunc = inds -> bathIntMatrix[inds...]
                    if numBathSites ≥ 2
                        cavityIndices = [numBathSites-1, numBathSites-3]
                    else
                        cavityIndices = Int64[]
                    end

                    hamiltonian = KondoModel(dispersion, kondoJMatrix, collect(1:numBathSites), bathIntFunc; bathIntLegs=bathIntLegs, globalField=globalField, cavityIndices=cavityIndices)
                    comparisonHamiltonian = []
                    for i in 0:numBathSites
                        push!(comparisonHamiltonian, ("n", [1 + 2 * i], globalField))
                        push!(comparisonHamiltonian, ("n", [2 + 2 * i], -globalField))
                    end
                    for i in 1:numBathSites
                        if abs(dispersion[i]) > 1e-15
                            push!(comparisonHamiltonian, ("n", [1 + 2 * i], dispersion[i]))
                            push!(comparisonHamiltonian, ("n", [2 + 2 * i], dispersion[i]))
                        end
                    end
                    for (i,j) in Iterators.product(1:numBathSites, 1:numBathSites)
                        if i ∈ cavityIndices || j ∈ cavityIndices
                            continue
                        end
                        if abs(kondoJMatrix[i, j]) > 1e-15
                            push!(comparisonHamiltonian, ("+-+-", [1, 2, 2 + 2 * i, 1 + 2 * j], kondoJMatrix[i, j]/2))
                            push!(comparisonHamiltonian, ("+-+-", [2, 1, 1 + 2 * i, 2 + 2 * j], kondoJMatrix[i, j]/2))
                            push!(comparisonHamiltonian, ("n+-", [1, 1 + 2 * i, 1 + 2 * j], kondoJMatrix[i, j]/4))
                            push!(comparisonHamiltonian, ("n+-", [1, 2 + 2 * i, 2 + 2 * j], -kondoJMatrix[i, j]/4))
                            push!(comparisonHamiltonian, ("n+-", [2, 1 + 2 * i, 1 + 2 * j], -kondoJMatrix[i, j]/4))
                            push!(comparisonHamiltonian, ("n+-", [2, 2 + 2 * i, 2 + 2 * j], kondoJMatrix[i, j]/4))
                        end
                    end

                    for (i,j, k, l) in Iterators.product(repeat([1:numBathSites], 4)...)
                        bathInt = bathIntMatrix[i, j, k, l]
                        if length(unique([i, j, k, l])) > bathIntLegs || abs(bathInt) < 1e-15
                            continue
                        end
                        push!(comparisonHamiltonian, ("+-+-", [2 * i + 1, 2 * j + 1, 2 * k + 1, 2 * l + 1], -bathInt/2))
                        push!(comparisonHamiltonian, ("+-+-", [2 * i + 2, 2 * j + 2, 2 * k + 2, 2 * l + 2], -bathInt/2))
                        push!(comparisonHamiltonian, ("+-+-", [2 * i + 1, 2 * j + 1, 2 * k + 2, 2 * l + 2], bathInt/2))
                        push!(comparisonHamiltonian, ("+-+-", [2 * i + 2, 2 * j + 2, 2 * k + 1, 2 * l + 1], bathInt/2))
                    end

                    @test issetequal(hamiltonian, comparisonHamiltonian)
                end
            end
        end
    end
end
