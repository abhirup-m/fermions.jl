@testset "Groundstate Correlation" begin
    basisStates = BasisStates(4)
    hop_t = rand()
    U = rand()
    eps = -U/2
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


@testset "Entanglement entropy" begin
    @testset "θ = $(theta)" for theta in rand(10) .* 2π
        state = Dict(BitVector([1, 0]) => cos(theta), BitVector([0, 1]) => sin(theta))
        @test vnEntropy(state, [1]) ≈ vnEntropy(state, [2]) ≈ -cos(theta)^2 * log(cos(theta)^2) - sin(theta)^2 * log(sin(theta)^2)
    end

    coeffs = ones(3)
    state = Dict(BitVector([1, 0, 1]) => coeffs[1], BitVector([1, 1, 0]) => coeffs[2], BitVector([0, 1, 1]) => coeffs[3])
    @test vnEntropy(state, [3]) ≈ log(3)/3 + log(3/2) * 2 / 3
end
