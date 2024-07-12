@testset "Groundstate Correlation" begin
    basisStates = BasisStates(4)
    hop_t = rand()
    U = rand()
    eps = -U / 2
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


@testset "Entanglement entropy, Mutual information" begin

    @testset "θ = $(theta)" for theta in rand(10) .* 2π
        state = Dict(BitVector([1, 0]) => cos(theta), BitVector([0, 1]) => sin(theta))
        @test vnEntropy(state, [1]) ≈ vnEntropy(state, [2]) ≈ -cos(theta)^2 * log(cos(theta)^2) - sin(theta)^2 * log(sin(theta)^2)
    end

    coeffs = ones(3)
    state = Dict(BitVector([1, 0, 1]) => coeffs[1], BitVector([1, 1, 0]) => coeffs[2], BitVector([0, 1, 1]) => coeffs[3])
    @test vnEntropy(state, [1, 2]) ≈ vnEntropy(state, [3]) ≈ log(3) / 3 + log(3 / 2) * 2 / 3

    state = Dict(BitVector([1, 0, 1]) => rand(), BitVector([0, 1, 1]) => rand())
    @test vnEntropy(state, [3]) ≈ 0

    coeffs = [1.0, 1.0, 1.0]#[rand() for _ in 1:3]
    state = Dict(BitVector([1, 0, 1]) => coeffs[1], BitVector([0, 1, 1]) => coeffs[2], BitVector([1, 1, 0]) => coeffs[3])
    coeffs ./= sum(coeffs .^ 2)^0.5
    @test vnEntropy(state, [1]) ≈ -(coeffs[1]^2 + coeffs[3]^2) * log(coeffs[1]^2 + coeffs[3]^2) - (coeffs[2]^2) * log(coeffs[2]^2)
    @test vnEntropy(state, [2]) ≈ -(coeffs[3]^2 + coeffs[2]^2) * log(coeffs[3]^2 + coeffs[2]^2) - (coeffs[1]^2) * log(coeffs[1]^2)
    @test vnEntropy(state, [1, 2]) ≈ -(coeffs[1]^2 + coeffs[2]^2) * log(coeffs[1]^2 + coeffs[2]^2) - (coeffs[3]^2) * log(coeffs[3]^2)
    @test mutInfo(state, ([1], [2])) == vnEntropy(state, [1]) + vnEntropy(state, [2]) - vnEntropy(state, [1, 2])

end

@testset "Spectral Function" begin
    basisStates = BasisStates(4)
    hop_t = abs(rand())
    U = abs(rand())
    eps = -U / 2
    broadening = 1e-3
    operatorList = HubbardDimerOplist([eps, eps, eps, eps], [U, U], [hop_t, hop_t])
    eigvals, eigvecs = Spectrum(operatorList, basisStates)
    omegaVals = collect(range(-10.0, stop=10.0, length=1000))
    specfunc = SpecFunc(eigvals, eigvecs, ("-", [1], 1.0), ("+", [1], 1.0), omegaVals, broadening)
    specfuncCompare = HubbardDimerSpecFunc(eps, U, hop_t, omegaVals, broadening)
    @test specfunc ./ maximum(specfunc) ≈ specfuncCompare ./ maximum(specfuncCompare)
end
