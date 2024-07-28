@testset "Spectrum" begin
    eps = rand()
    hop_t = rand()
    U = rand()
    @testset "sector=$((n, m))" for (n, m) in [(0, 0), (1, 1), (1, -1), (2, 2), (2, 0), (2, -2), (3, 1), (3, -1), (4, 0)]
        basis = BasisStates(4; totOccReq=n, magzReq=m)
        operatorList = HubbardDimerOplist(eps, U, hop_t)
        eigvals, eigvecs = Spectrum(operatorList, basis)
        comparisonMatrix = HubbardDimerMatrix(eps, U, hop_t)[(n, m)]
        eigvalTest, eigvecTest = eigen(comparisonMatrix)
        @test eigvals ≈ eigvalTest
        for (i, v2) in enumerate(eachcol(eigvecTest))
            v2 = v2[abs.(v2) .> 1e-14]
            if (n,m) == (2, 0) && BitVector([1, 0, 0, 1]) in keys(eigvecs[i]) && BitVector([1, 1, 0, 0]) in keys(eigvecs[i])
                vals = [eigvecs[i][BitVector([1, 0, 0, 1])], eigvecs[i][BitVector([0, 1, 1, 0])], eigvecs[i][BitVector([1, 1, 0, 0])], eigvecs[i][BitVector([0, 0, 1, 1])]]
                vals = vals[abs.(vals) .> 1e-14]
                @test vals ./ maximum(abs.(vals)) ≈ v2 ./ maximum(abs.(v2)) || vals ./ maximum(abs.(vals)) ≈ -1 .* v2 ./ maximum(abs.(v2))
            else
                @test collect(values(eigvecs[i])) ./ collect(values(eigvecs[i]))[1] ≈ v2 ./ v2[1]
            end
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
    operatorList = HubbardDimerOplist(eps, U, hop_t)
    eigvals, eigvecs = Spectrum(operatorList, basisStates)
    @test eigvals[1] ≈ 2 * eps + U / 2 - Δ / 2
    @test Set(keys(eigvecs[1])) == Set([[1, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 0]])
    scaledGstate = Dict(k => v / eigvecs[1][[1, 0, 0, 1]] for (k, v) in eigvecs[1])
    @test scaledGstate[[1, 0, 0, 1]] ≈ 1.0
    @test scaledGstate[[0, 1, 1, 0]] ≈ -1.0
    @test scaledGstate[[0, 0, 1, 1]] ≈ a1 / a2
    @test scaledGstate[[1, 1, 0, 0]] ≈ a1 / a2

end


@testset "Degenerate Ground States" begin
    basisStates = BasisStates(4)
    eps = -abs(rand())
    hop_t = 0
    U = -2 * eps
    operatorList = HubbardDimerOplist(eps, U, hop_t)
    eigvals, eigvecs = Spectrum(operatorList, basisStates)
    @test eigvals[1] == eigvals[2] == eigvals[3] == eigvals[4] ≈ 2 * eps
    @test Set((eigvecs[1], eigvecs[2], eigvecs[3], eigvecs[4])) == Set((Dict([1, 0, 1, 0] => 1.0), Dict([0, 1, 0, 1] => 1.0), Dict([0, 1, 1, 0] => 1.0), Dict([1, 0, 0, 1] => 1.0)))
end
