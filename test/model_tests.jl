@testset "J-Ub model, matrix" begin
    Ek, J, W = rand(3)
    dispDict = Dict(1 => Ek)
    kondoDict = Dict((1,1) => J)
    bathIntDict = Dict((1,1,1,1) => W)
    basis = BasisStates(4)
    oplist = kondoKSpace(dispDict, kondoDict, bathIntDict)
    matrixBlocks = generalOperatorMatrix(basis, oplist)
    @test matrixBlocks[(0,0)] == [0.0;;]
    @test matrixBlocks[(4,0)] == [2 * Ek;;]
    @test matrixBlocks[(1,1)] == matrixBlocks[(1,-1)] == [[Ek - W/2, 0] [0, 0]]
    @test matrixBlocks[(3,1)] == matrixBlocks[(3,1)] == [[2*Ek, 0] [0, Ek - W/2]]
    @test matrixBlocks[(2, 0)] == [[2*Ek, 0, 0, 0] [0, Ek - J/4 - W/2, J/2, 0] [0, J/2, Ek - J/4 - W/2, 0] [0, 0, 0, 0]]
    @test matrixBlocks[(2, 2)] ≈ [Ek + J/4 - W/2;;]
    @test matrixBlocks[(2, -2)] ≈ [Ek + J/4 - W/2;;]
end

@testset "J-Ub model, symmetry" begin
    num = 3
    basis = BasisStates(2 * num + 2)
    occupancyOp = i -> Dict{Tuple{String,Vector{Int64}},Float64}(("n", [2 * i + 1]) => 1.0)
    spinFlipOp = i -> Dict{Tuple{String,Vector{Int64}},Float64}(("+-+-", [1, 2, 2 * i + 2, 2 * i + 1]) => 1.0)
    Ek = zeros(num)
    dispDict = Dict(zip(1:num, Ek))
    kvals = [(0, -π), (π/2, -π/2), (π, 0)]
    kondoDict = Dict((i,j) => cos(kvals[i][1] - kvals[j][1]) + cos(kvals[i][2] - kvals[j][2]) for (i,j) in Iterators.product(1:num, 1:num))
    bathIntDict = Dict((i,j,k,l) => cos(kvals[i][1] - kvals[j][1] + kvals[k][1] - kvals[l][1]) + cos(kvals[i][2] - kvals[j][2] + kvals[k][2] - kvals[l][2]) 
                     for (i,j,k,l) in Iterators.product(1:num, 1:num, 1:num, 1:num))
    eigvals, eigvecs = getSpectrum(generalOperatorMatrix(basis, kondoKSpace(dispDict, kondoDict, bathIntDict)))
    correlationOcc = gstateCorrelation(basis, eigvals, eigvecs, occupancyOp.(1:num))
    correlationSpinFlip = gstateCorrelation(basis, eigvals, eigvecs, spinFlipOp.(1:num))
    @test correlationOcc[1] ≈ correlationOcc[3]
    @test correlationSpinFlip[1] ≈ correlationSpinFlip[3]
    kvals = reverse([(0, -π), (π/2, -π/2), (π, 0)])
    kondoDict = Dict((i,j) => cos(kvals[i][1] - kvals[j][1]) + cos(kvals[i][2] - kvals[j][2]) for (i,j) in Iterators.product(1:num, 1:num))
    bathIntDict = Dict((i,j,k,l) => cos(kvals[i][1] - kvals[j][1] + kvals[k][1] - kvals[l][1]) + cos(kvals[i][2] - kvals[j][2] + kvals[k][2] - kvals[l][2]) 
                     for (i,j,k,l) in Iterators.product(1:num, 1:num, 1:num, 1:num))
    eigvals, eigvecs = getSpectrum(generalOperatorMatrix(basis, kondoKSpace(dispDict, kondoDict, bathIntDict)))
    revcorrelationOcc = gstateCorrelation(basis, eigvals, eigvecs, occupancyOp.(1:num))
    revcorrelationSpinFlip = gstateCorrelation(basis, eigvals, eigvecs, spinFlipOp.(1:num))
    @test revcorrelationOcc == correlationOcc
    @test revcorrelationSpinFlip == correlationSpinFlip
end
