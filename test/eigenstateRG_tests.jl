@testset "Eigenstate RG" begin
    initState = Dict(BitVector([1, 0, 0, 1]) => 1 / √2, BitVector([0, 1, 1, 0]) => -1 / √2)
    numSteps = 8
    alphaValues = ones(numSteps)
    unitaryOperatorFunction = simpleUnitaries
    stateExpansionFunction = stateExpFunc1Ck
    sectors = "ph"
    @time stateFlow = getWavefunctionRG(initState, alphaValues, numSteps, unitaryOperatorFunction, stateExpansionFunction, sectors)
    @assert false
    @test stateFlow[1] == initState
    scaledState = Dict(k => v / stateFlow[2][[1, 0, 0, 1, 1, 1, 0, 0]] for (k, v) in stateFlow[2])
    @test Set(keys(scaledState)) == Set([BitVector([1, 0, 0, 1, 1, 1, 0, 0]), BitVector([1, 0, 1, 1, 0, 1, 0, 0]), BitVector([1, 0, 0, 0, 1, 1, 0, 1]), BitVector([0, 1, 1, 0, 1, 1, 0, 0]), BitVector([0, 1, 1, 1, 1, 0, 0, 0]), BitVector([0, 1, 0, 0, 1, 1, 1, 0])])
    @test scaledState[BitVector([1, 0, 0, 1, 1, 1, 0, 0])] ≈ 1.0
    @test scaledState[BitVector([1, 0, 1, 1, 0, 1, 0, 0])] ≈ -3 / 4
    @test scaledState[BitVector([1, 0, 0, 0, 1, 1, 0, 1])] ≈ -3 / 4
    @test scaledState[BitVector([0, 1, 1, 0, 1, 1, 0, 0])] ≈ -1.0
    @test scaledState[BitVector([0, 1, 1, 1, 1, 0, 0, 0])] ≈ 3 / 4
    @test scaledState[BitVector([0, 1, 0, 0, 1, 1, 1, 0])] ≈ 3 / 4
end

