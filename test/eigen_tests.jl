@testset "Spectrum" begin
    eps = rand(4)
    hop_t = rand(2)
    U = rand(2)
    @testset "sector=$((n, m))" for (n, m) in [(0, 0), (1, 1), (1, -1), (2, 2), (2, 0), (2, -2), (3, 1), (3, -1), (4, 0)]
        basis = BasisStates(4; totOccCriteria=(x, N) -> sum(x) == n, magzCriteria=x -> sum(x[1:2:end]) - sum(x[2:2:end]) == m)
        operatorList = HubbardDimerOplist(eps, U, hop_t)
        eigvals, eigvecs = Spectrum(operatorList, basis)
        comparisonMatrix = HubbardDimerMatrix(eps, U, hop_t)[(n, m)]
        eigvalTest, eigvecTest = eigen(comparisonMatrix)
        if (n, m) == (2, 0)
            display(eigvecs)
            display(eigvecTest)
        end
        @test eigvals ≈ eigvalTest
        for (v1, v2) in zip(eigvecs, [eigvecTest[:, i] for i in eachindex(eachrow(eigvecTest))])
            @test collect(values(v1)) ≈ collect(v2)
        end
    end
end


@testset "Ground States" begin
    basisStates = BasisStates(4)
    eps = rand()
    hop_t = rand()
    U = rand()
    Δ = (U^2 + 16 * hop_t^2)^0.5
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    eigvals, eigvecs = Spectrum(operatorList, basisStates)
    @test eigvals[1] ≈ 2 * eps + U / 2 - Δ / 2
end
@assert false

@testset "Degenerate Ground States" begin
    basisStates = BasisStates(4)
    eps = -abs(rand())
    hop_t = 0
    U = -2 * eps
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    computedMatrix = generalOperatorMatrix(basisStates, operatorList)
    eigvals, eigvecs = Spectrum(computedMatrix)
    gsEnergy, groundStates, blocks = getGstate(eigvals, eigvecs)
    @test gsEnergy ≈ 2 * eps
    @test Set(blocks) == Set([(2, 0), (2, 0), (2, 2), (2, -2)])
    @test sort(join.(groundStates)) == sort(join.([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [1.0], [1.0]]))
end
