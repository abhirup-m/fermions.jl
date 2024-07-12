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
        @test eigvals ≈ eigvalTest
        for (v1, v2) in zip(eigvecs, [eigvecTest[:, i] for i in eachindex(eachrow(eigvecTest))])
            @test collect(values(v1)) ≈ collect(v2)
        end
    end
end


@testset "Ground States" begin
    basisStates = BasisStates(4)
    hop_t = abs(rand())
    U = abs(rand())
    eps = -U / 2
    Δ = (U^2 + 16 * hop_t^2)^0.5
    a1 = 0.5 * √((Δ - U) / Δ)
    a2 = 2 * hop_t / √(Δ * (Δ - U))
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    eigvals, eigvecs = Spectrum(operatorList, basisStates)
    @test eigvals[1] ≈ 2 * eps + U / 2 - Δ / 2
    @test Set(keys(eigvecs[1])) == Set([[1, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 0]])
    scaledGstate = Dict(k => v / eigvecs[1][[1, 0, 0, 1]] for (k, v) in eigvecs[1])
    @test scaledGstate[[1, 0, 0, 1]] == 1.0
    @test scaledGstate[[0, 1, 1, 0]] == -1.0
    @test scaledGstate[[0, 0, 1, 1]] ≈ a1 / a2
    @test scaledGstate[[1, 1, 0, 0]] ≈ a1 / a2

end


@testset "Degenerate Ground States" begin
    basisStates = BasisStates(4)
    eps = -abs(rand())
    hop_t = 0
    U = -2 * eps
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    eigvals, eigvecs = Spectrum(operatorList, basisStates)
    @test eigvals[1] == eigvals[2] == eigvals[3] == eigvals[4] ≈ 2 * eps
    @test Set((eigvecs[1], eigvecs[2], eigvecs[3], eigvecs[4])) == Set((Dict([1, 0, 1, 0] => 1.0), Dict([0, 1, 0, 1] => 1.0), Dict([0, 1, 1, 0] => 1.0), Dict([1, 0, 0, 1] => 1.0)))
end
