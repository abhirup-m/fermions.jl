@testset "Groundstate Correlation" begin
    basisStates = BasisStates(4)
    eps = rand()
    hop_t = rand()
    U = rand()
    Δ = (U^2 + 16*hop_t^2)^0.5
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    computedMatrix = generalOperatorMatrix(basisStates, operatorList)
    eigvals, eigvecs = getSpectrum(computedMatrix)
    doubOccOperator = Dict(("nn", [1,2]) => 0.5, ("nn", [3,4]) => 0.5)
    totSzOperator = Dict(("n", [1]) => 0.5, ("n", [2]) => -0.5, ("n", [3]) => 0.5, ("n", [4]) => -0.5)
    spinFlipOperator = Dict(("+-+-", [1,2,4,3]) => 1.0, ("+-+-", [2,1,3,4]) => 1.0)
    doubOccValue = gstateCorrelation(basisStates, eigvals, eigvecs, doubOccOperator)
    totSzValue = gstateCorrelation(basisStates, eigvals, eigvecs, totSzOperator)
    spinFlipValue = gstateCorrelation(basisStates, eigvals, eigvecs, spinFlipOperator)
    @test doubOccValue ≈ (Δ - U) / (4 * Δ)
    @test totSzValue ≈ 0
    @test spinFlipValue ≈ 8 * hop_t^2 / (Δ * (Δ - U))
end
