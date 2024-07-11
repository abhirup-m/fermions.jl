@testset "Groundstate Correlation" begin
    basisStates = BasisStates(4)
    eps = rand()
    hop_t = rand()
    U = rand()
    Δ = (U^2 + 16 * hop_t^2)^0.5
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    eigvals, eigvecs = Spectrum(operatorList, basisStates)

    doubOccOperator = [("nn", [1, 2], 0.5), ("nn", [3, 4], 0.5)]
    totSzOperator = [("n", [1], 0.5), ("n", [2], -0.5), ("n", [3], 0.5), ("n", [4], -0.5)]
    spinFlipOperator = [("+-+-", [1, 2, 4, 3], 1.0), ("+-+-", [3, 4, 2, 1], 1.0)]

    @test FastCorrelation(eigvecs[1], doubOccOperator) ≈ (Δ - U) / (4 * Δ)
    @test FastCorrelation(eigvecs[1], totSzOperator) ≈ 0
    @test FastCorrelation(eigvecs[1], spinFlipOperator) ≈ -8 * hop_t^2 / (Δ * (Δ - U))
end
