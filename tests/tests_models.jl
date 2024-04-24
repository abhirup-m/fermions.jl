using Test
include("../source/models.jl")

@testset "Kondo k-space model" begin
    function bathIntFunc(points)
        kx = [[-pi, -pi, 0, 0][p] for p in points]
        ky = [[-pi, 0, -pi, 0][p] for p in points]
        return 0.5 .* (cos(kx[1] - kx[2] + kx[3] - kx[4]) + cos(ky[1] - ky[2] + ky[3] - ky[4]))
    end
    chosenIndices = [1, 2]
    size_BZ = 2
    basis = BasisStates(length(chosenIndices) * 2 + 2)
    kondoJArray = Array{Float64}(undef, size_BZ^2, size_BZ^2)
    for (p1, p2) in Iterators.product(1:size_BZ^2, 1:size_BZ^2)
        k1x, k1y = [(-pi, -pi), (-pi, 0), (0, -pi), (0, 0)][p1]
        k2x, k2y = [(-pi, -pi), (-pi, 0), (0, -pi), (0, 0)][p2]
        kondoJArray[p1, p2, 1] = 0.5 .* (cos.(k1x .- k2x) .+ cos.(k1y .- k2y))
    end

    dispersionArray = -2.0 .* [-2, 0, 0, 2]
    oplist = KondoKSpace(chosenIndices, dispersionArray, kondoJArray, bathIntFunc)
    matrices = generalOperatorMatrix(basis, oplist)
    @testset "(n, m) = $key" for key in keys(matrices)
        @test ishermitian(matrices[key])
    end
end

