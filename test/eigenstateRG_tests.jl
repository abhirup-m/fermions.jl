@testset "Eigenstate RG 1CK ph N=1" begin
    initState = Dict(BitVector([1, 0, 0, 1]) => 1 / √2, BitVector([0, 1, 1, 0]) => -1 / √2)
    numSteps = 1
    alphaValues = ones(numSteps)
    sectors = "ph"
    stateFlow = getWavefunctionRG(initState, alphaValues, numSteps, unitaries1CK, stateExpansion1CK, sectors)

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

@testset "Eigenstate RG 1CK ph N=2" begin
    initState = Dict(
        BitVector([1, 0, 0, 1, 0, 0]) => 1.0,
        BitVector([1, 0, 0, 0, 0, 1]) => 1.0,
        BitVector([0, 1, 1, 0, 0, 0]) => -1.0,
        BitVector([0, 1, 0, 0, 1, 0]) => -1.0,
    )
    refState = [1, 0, 0, 1, 0, 0, 1, 1, 0, 0]
    numSteps = 1
    alphaValues = ones(numSteps)
    sectors = "ph"
    stateFlow = getWavefunctionRG(initState, alphaValues, numSteps, unitaries1CK, stateExpansion1CK, sectors)
    @test stateFlow[1] == initState

    scaledState = Dict(k => v / stateFlow[2][refState] for (k, v) in stateFlow[2])
    @test Set(keys(scaledState)) == Set([
        [1, 0, 0, 1, 0, 0, 1, 1, 0, 0],
        [0, 1, 1, 0, 0, 0, 1, 1, 0, 0],
        [1, 0, 0, 0, 1, 1, 0, 1, 0, 0],
        [0, 1, 0, 0, 1, 1, 1, 0, 0, 0],
        [1, 0, 0, 0, 0, 1, 1, 1, 0, 0],
        [0, 1, 0, 0, 1, 0, 1, 1, 0, 0],
        [1, 0, 1, 1, 0, 0, 0, 1, 0, 0],
        [0, 1, 1, 1, 0, 0, 1, 0, 0, 0],
        [1, 0, 1, 0, 0, 1, 0, 1, 0, 0],
        [0, 1, 0, 1, 1, 0, 1, 0, 0, 0],
        [1, 0, 0, 1, 1, 0, 0, 1, 0, 0],
        [0, 1, 1, 0, 0, 1, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 1, 1, 1, 0],
        [1, 0, 0, 0, 0, 0, 1, 1, 0, 1],
    ])

    @test scaledState[[1, 0, 0, 1, 0, 0, 1, 1, 0, 0]] == 1.0
    @test scaledState[[0, 1, 1, 0, 0, 0, 1, 1, 0, 0]] == -1.0
    @test scaledState[[1, 0, 0, 0, 1, 1, 0, 1, 0, 0]] == -0.75
    @test scaledState[[0, 1, 0, 0, 1, 1, 1, 0, 0, 0]] == 0.75
    @test scaledState[[1, 0, 0, 0, 0, 1, 1, 1, 0, 0]] == 1.0
    @test scaledState[[0, 1, 0, 0, 1, 0, 1, 1, 0, 0]] == -1.0
    @test scaledState[[1, 0, 1, 1, 0, 0, 0, 1, 0, 0]] == -0.75
    @test scaledState[[0, 1, 1, 1, 0, 0, 1, 0, 0, 0]] == 0.75
    @test scaledState[[1, 0, 1, 0, 0, 1, 0, 1, 0, 0]] == -0.75
    @test scaledState[[0, 1, 0, 1, 1, 0, 1, 0, 0, 0]] == -0.75
    @test scaledState[[1, 0, 0, 1, 1, 0, 0, 1, 0, 0]] == 0.75
    @test scaledState[[0, 1, 1, 0, 0, 1, 1, 0, 0, 0]] == 0.75
    @test scaledState[[0, 1, 0, 0, 0, 0, 1, 1, 1, 0]] == 1.5
    @test scaledState[[1, 0, 0, 0, 0, 0, 1, 1, 0, 1]] == -1.5
end

@testset "Eigenstate RG 1CK p N=2" begin
    initState = Dict(
        BitVector([1, 0, 0, 1, 0, 0]) => 1.0,
        BitVector([1, 0, 0, 0, 0, 1]) => 1.0,
        BitVector([0, 1, 1, 0, 0, 0]) => -1.0,
        BitVector([0, 1, 0, 0, 1, 0]) => -1.0,
    )
    refState = [1, 0, 0, 1, 0, 0, 1, 1]
    numSteps = 1
    alphaValues = ones(numSteps)
    sectors = "p"
    stateFlow = getWavefunctionRG(initState, alphaValues, numSteps, unitaries1CK, stateExpansion1CK, sectors)
    @test stateFlow[1] == initState

    scaledState = Dict(k => v / stateFlow[2][refState] for (k, v) in stateFlow[2])
    @test Set(keys(scaledState)) == Set([
        [1, 0, 0, 1, 0, 0, 1, 1],
        [0, 1, 1, 0, 0, 0, 1, 1],
        [1, 0, 0, 0, 1, 1, 0, 1],
        [0, 1, 0, 0, 1, 1, 1, 0],
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [1, 0, 1, 1, 0, 0, 0, 1],
        [0, 1, 1, 1, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 1, 0, 1],
        [0, 1, 0, 1, 1, 0, 1, 0],
        [1, 0, 0, 1, 1, 0, 0, 1],
        [0, 1, 1, 0, 0, 1, 1, 0],
    ])

    @test scaledState[[1, 0, 0, 1, 0, 0, 1, 1]] == 1.0
    @test scaledState[[0, 1, 1, 0, 0, 0, 1, 1]] == -1.0
    @test scaledState[[1, 0, 0, 0, 1, 1, 0, 1]] == -0.75
    @test scaledState[[0, 1, 0, 0, 1, 1, 1, 0]] == 0.75
    @test scaledState[[1, 0, 0, 0, 0, 1, 1, 1]] == 1.0
    @test scaledState[[0, 1, 0, 0, 1, 0, 1, 1]] == -1.0
    @test scaledState[[1, 0, 1, 1, 0, 0, 0, 1]] == -0.75
    @test scaledState[[0, 1, 1, 1, 0, 0, 1, 0]] == 0.75
    @test scaledState[[1, 0, 1, 0, 0, 1, 0, 1]] == -0.75
    @test scaledState[[0, 1, 0, 1, 1, 0, 1, 0]] == -0.75
    @test scaledState[[1, 0, 0, 1, 1, 0, 0, 1]] == 0.75
    @test scaledState[[0, 1, 1, 0, 0, 1, 1, 0]] == 0.75
end

@testset "Eigenstate RG 1CK h N=2" begin
    initState = Dict(
        BitVector([1, 0, 0, 1, 0, 0]) => 1.0,
        BitVector([1, 0, 0, 0, 0, 1]) => 1.0,
        BitVector([0, 1, 1, 0, 0, 0]) => -1.0,
        BitVector([0, 1, 0, 0, 1, 0]) => -1.0,
    )
    refState = [1, 0, 0, 1, 0, 0, 0, 0]
    numSteps = 1
    alphaValues = ones(numSteps)
    sectors = "h"
    stateFlow = getWavefunctionRG(initState, alphaValues, numSteps, unitaries1CK, stateExpansion1CK, sectors)
    @test stateFlow[1] == initState

    scaledState = Dict(k => v / stateFlow[2][refState] for (k, v) in stateFlow[2])
    @test Set(keys(scaledState)) == Set([
        [1, 0, 0, 1, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0, 0, 0, 1],
    ])

    @test scaledState[[1, 0, 0, 1, 0, 0, 0, 0]] == 1.0
    @test scaledState[[0, 1, 1, 0, 0, 0, 0, 0]] == -1.0
    @test scaledState[[1, 0, 0, 0, 0, 1, 0, 0]] == 1.0
    @test scaledState[[0, 1, 0, 0, 1, 0, 0, 0]] == -1.0
    @test scaledState[[0, 1, 0, 0, 0, 0, 1, 0]] == 1.5
    @test scaledState[[1, 0, 0, 0, 0, 0, 0, 1]] == -1.5
end

@testset "Eigenstate RG linearity" begin
    bstates = [
        BitVector([1, 0, 0, 1]),
        BitVector([0, 1, 1, 0]),
    ]
    initState = Dict(bstates .=> 1.0)
    numSteps = 1
    alphaValues = rand(numSteps)
    sectors = "ph"
    totalFinalState = getWavefunctionRG(initState, alphaValues, numSteps, unitaries1CK, stateExpansion1CK, sectors)[end]
    scaledTotalFinal = Dict(keys(totalFinalState) .=> values(totalFinalState) ./ maximum(values(totalFinalState)))

    concoctedFinalState = mergewith(+, [getWavefunctionRG(Dict(k => 1.0), alphaValues, numSteps, unitaries1CK, stateExpansion1CK, sectors)[end]
                                        for (k, _) in initState]...)
    scaledConcoctedFinal = Dict(k => v / maximum(values(concoctedFinalState)) for (k, v) in concoctedFinalState)
    @test Set(keys(scaledTotalFinal)) == Set(keys(scaledConcoctedFinal))
    @test all(key -> scaledTotalFinal[key] ≈ scaledConcoctedFinal[key], keys(scaledTotalFinal))
end
