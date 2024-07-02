@testset "Expand basis" begin
    t = 1.0
    U = 0.0
    initBasis = BasisStates(4)
    retainSize = 100
    hamiltonianFamily = [dimerHamiltonian(U, t), dimerAdditionalHamiltonian(U, t, 3)]
    numStatesFamily = Int64[2, 3]
    hamiltonianMatrix = operatorMatrix(initBasis, hamiltonianFamily[1])

    spectrum = getSpectrum(initBasis, hamiltonianMatrix; maxNum=retainSize)

    basisStates = spectrum[2]
    basisStates, diagElements = expandBasis(basisStates, numStatesFamily[2] - numStatesFamily[1], spectrum[1], retainSize)

    @test isempty(setdiff(keys(basisStates), [(0, 0), (1/6, -1), (1/6, 1), (2/6, 0), (4/6, 0), (5/6, 1), (5/6, -1), (4/6, 2), (4/6, -2), (6/6, 0), (3/6, 3), (3/6, -3), (2/6, 2), (2/6, -2), (3/6, 1), (3/6, -1), (4/6, 2), (4/6, -2)]))
    @test basisStates[(0, 0)] == [Dict([0, 0, 0, 0, 0, 0] => 1.0)]
    @test basisStates[(1/6, -1)] == [Dict([0, 0, 0, 0, 0, 1] => 1.0), 
                                     Dict([0, 1, 0, 0, 0, 0] => -1/2^0.5, [0, 0, 0, 1, 0, 0] => -1/2^0.5), 
                                     Dict([0, 1, 0, 0, 0, 0] => -1/2^0.5, [0, 0, 0, 1, 0, 0] => 1/2^0.5)
                                    ]
    # @test basisStates[(1/6, 1)] == 
    # @test basisStates[(2/6, 0)] == 
    # @test basisStates[(4/6, 0)] == 
    # @test basisStates[(5/6, 1)] == 
    # @test basisStates[(5/6, -1)] == 
    # @test basisStates[(4/6, 2)] == 
    # @test basisStates[(4/6, -2)] == 
    # @test basisStates[(6/6, 0)] == 
    # @test basisStates[(3/6, 3)] == 
    # @test basisStates[(3/6, -3)] == 
    # @test basisStates[(2/6, 2)] == 
    # @test basisStates[(2/6, -2)] == 
    # @test basisStates[(3/6, 1)] == 
    # @test basisStates[(3/6, -1)] == 
    # @test basisStates[(4/6, 2)] == 
    # @test basisStates[(4/6, -2)] == 
end
