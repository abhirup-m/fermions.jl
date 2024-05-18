@testset "Spectrum" begin
    basisStates = BasisStates(4)
    eps = rand(4)
    hop_t = rand(2)
    U = rand(2)
    operatorList = HubbardDimerOplist(eps, U, hop_t)
    computedMatrix = generalOperatorMatrix(basisStates, operatorList)
    eigvals, eigvecs = getSpectrum(computedMatrix)
    comparisonMatrix = HubbardDimerMatrix(eps, U, hop_t)
    @test keys(comparisonMatrix) == keys(computedMatrix)
    for key in keys(comparisonMatrix)
        eigvalTest, eigvecTest = eigen(comparisonMatrix[key])
        @test eigvals[key] ≈ eigvalTest
        @test eigvecs[key] ≈ [eigvecTest[:, i] for i in eachindex(eachrow(eigvecTest))]
    end
end


@testset "Ground States" begin
    basisStates = BasisStates(4)
    eps = rand()
    hop_t = rand()
    U = rand()
    Δ = (U^2 + 16*hop_t^2)^0.5
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    computedMatrix = generalOperatorMatrix(basisStates, operatorList)
    eigvals, eigvecs = getSpectrum(computedMatrix)
    gsEnergy, groundStates = getGstate(eigvals, eigvecs)
    @test gsEnergy ≈ 2 * eps + U / 2 - Δ / 2
end

@testset "Degenerate Ground States" begin
    basisStates = BasisStates(4)
    eps = -abs(rand())
    hop_t = 0
    U = -2 * eps
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    computedMatrix = generalOperatorMatrix(basisStates, operatorList)
    eigvals, eigvecs = getSpectrum(computedMatrix)
    gsEnergy, groundStates, blocks = getGstate(eigvals, eigvecs)
    @test gsEnergy ≈ 2 * eps
    @test Set(blocks) == Set([(2, 0), (2,0), (2,2), (2,-2)])
    @test sort(join.(groundStates)) == sort(join.([[0.0,1.0,0.0,0.0], [0.0,0.0,1.0,0.0], [1.0], [1.0]]))
end
