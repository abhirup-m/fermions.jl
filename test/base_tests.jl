@testset "BasisStates" begin

    @test issetequal(BasisStates(2),
        [
            Dict([0, 0] => 1.0),
            Dict([1, 0] => 1.0),
            Dict([0, 1] => 1.0),
            Dict([1, 1] => 1.0),
        ]
    )
    @test issetequal(BasisStates(3),
        [
            Dict([0, 0, 0] => 1.0),
            Dict([1, 0, 0] => 1.0),
            Dict([0, 1, 0] => 1.0),
            Dict([1, 1, 0] => 1.0),
            Dict([0, 0, 1] => 1.0),
            Dict([1, 0, 1] => 1.0),
            Dict([0, 1, 1] => 1.0),
            Dict([1, 1, 1] => 1.0),
        ]
    )
end


@testset "roundTo" begin
    @test roundTo(1, 1e-10) == 1
    @test roundTo(1 + 1e-11, 1e-10) == 1
    @test roundTo(1e-11, 1e-10) == 0
    @test roundTo(1e-10, 1e-10) == 1e-10
    @test roundTo(-1e-10, 1e-10) == -1e-10
    @test roundTo(0, 1e-10) == 0
end



@testset "TransformBit" begin
    @test TransformBit(false, 'n') == (0, 0)
    @test TransformBit(false, 'h') == (0, 1)
    @test TransformBit(false, '+') == (1, 1)
    @test TransformBit(false, '-') == (0, 0)
    @test TransformBit(true, 'n') == (1, 1)
    @test TransformBit(true, 'h') == (1, 0)
    @test TransformBit(true, '+') == (1, 0)
    @test TransformBit(true, '-') == (0, 1)
end

@testset "ApplyOperator" begin
    # checking linearity
    basis = BasisStates(4)
    coeffs1 = rand(4^4)
    coeffs2 = rand(length(basis))
    allOperators = [[(o1 * o2 * o3 * o4, [1, 2, 3, 4], coeffs1[i])]
                    for (i, (o1, o2, o3, o4)) in enumerate(Iterators.product(repeat([["+", "-", "n", "h"]], 4)...))]
    allStates = [Dict(k => coeffs2[i] * v for (k,v) in dict) for (i, dict) in enumerate(basis)]
    totalOperator = vcat(allOperators...)
    totalState = mergewith(+, allStates...)
    @test ApplyOperator(totalOperator, totalState) == mergewith(+, [ApplyOperator(operator, totalState) for operator in allOperators]...)
    @test all(x -> x ≈ 0, [v1 - v2 for (v1, v2) in zip(values(ApplyOperator(totalOperator, totalState)), values(mergewith(+, [ApplyOperator(totalOperator, state) for state in allStates]...)))])
    end

    # quantitative checking through basis states
    state = Dict(BitVector([0, 0]) => 0.1)
    @testset "state = [0, 0], operator = $op" for op in ["+", "-", "n", "h"]
        oplist = Dict((op, [1]) => 0.5)
        if op == "+"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 0]) => 0.05)
        elseif op == "h"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 0]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        oplist = Dict((op, [2]) => 0.5)
        if op == "+"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 1]) => 0.05)
        elseif op == "h"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 0]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        @testset for op2 in ["+", "-", "n", "h"]
            oplist = Dict((op * op2, [1, 2]) => 0.5)
            if occursin("-", op * op2) || occursin("n", op * op2)
                @test isempty(applyOperatorOnState(state, oplist))
            elseif op * op2 == "++"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 1]) => 0.05)
            elseif op * op2 == "+h"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 0]) => 0.05)
            elseif op * op2 == "h+"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 1]) => 0.05)
            elseif op * op2 == "hh"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 0]) => 0.05)
            end
        end
    end

    state = Dict(BitVector([1, 1]) => 0.1)
    @testset "state = [1, 1], operator = $op" for op in ["+", "-", "n", "h"]
        oplist = Dict((op, [1]) => 0.5)
        if op == "-"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 1]) => 0.05)
        elseif op == "n"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 1]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        oplist = Dict((op, [2]) => 0.5)
        if op == "-"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 0]) => -0.05)
        elseif op == "n"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 1]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        @testset for op2 in ["+", "-", "n", "h"]
            oplist = Dict((op * op2, [1, 2]) => 0.5)
            if occursin("+", op * op2) || occursin("h", op * op2)
                @test isempty(applyOperatorOnState(state, oplist))
            elseif op * op2 == "--"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 0]) => -0.05)
            elseif op * op2 == "-n"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 1]) => 0.05)
            elseif op * op2 == "n-"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 0]) => -0.05)
            elseif op * op2 == "nn"
                @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 1]) => 0.05)
            end
        end
    end
end


"""
@testset "GeneralOperatorMatrix" begin
    basisStates = BasisStates(4)
    eps = rand(4)
    hop_t = rand(2)
    U = rand(2)
    operatorList = HubbardDimerOplist(eps, U, hop_t)
    computedMatrix = generalOperatorMatrix(basisStates, operatorList)
    comparisonMatrix = HubbardDimerMatrix(eps, U, hop_t)
    @test keys(comparisonMatrix) == keys(computedMatrix)
    for key in keys(comparisonMatrix)
        @test computedMatrix[key] ≈ comparisonMatrix[key]
    end
end


@testset "GeneralOperatorMatrixSet" begin
    basisStates = BasisStates(4)
    eps_arr = [rand(4) for _ in 1:3]
    hop_t_arr = [rand(2) for _ in 1:3]
    U_arr = [rand(2) for _ in 1:3]
    couplingMatrix = Vector{Float64}[]
    operatorList = nothing
    for (eps, hop_t, U) in zip(eps_arr, hop_t_arr, U_arr)
        operatorTermsValues = HubbardDimerOplist(eps, U, hop_t)
        operatorList = collect(keys(operatorTermsValues))
        push!(couplingMatrix, collect(values(operatorTermsValues)))
    end
    computedMatrixSet = generalOperatorMatrix(basisStates, operatorList, couplingMatrix)
    comparisonMatrixSet = [HubbardDimerMatrix(eps_arr[i], U_arr[i], hop_t_arr[i]) for i in 1:3]
    @test length(computedMatrixSet) == length(comparisonMatrixSet)
    for (comparisonMatrix, computedMatrix) in zip(comparisonMatrixSet, computedMatrixSet)
        @test keys(comparisonMatrix) == keys(computedMatrix)
        for key in keys(comparisonMatrix)
            @test computedMatrix[key] ≈ comparisonMatrix[key]
        end
    end
end

@testset "transformBasis" begin
    basisStates = BasisStates(4)
    transformation = Dict((1, 1) => [[0.5, 0.5], [0.5, -0.5]], (2, 0) => [[0, 0, 2^0.5, 2^0.5], [0, 0, 2^0.5, -2^0.5], [1, 0, 0, 0], [0, 1, 0, 0]])
    transformedBasis = transformBasis(basisStates, transformation)
    for (k, v) in transformedBasis
        if k == (1, 1)
            @test transformedBasis[k][1] == Dict([1, 0, 0, 0] => 0.5, [0, 0, 1, 0] => 0.5)
            @test transformedBasis[k][2] == Dict([1, 0, 0, 0] => -0.5, [0, 0, 1, 0] => 0.5)
        elseif k == (2, 0)
            @test transformedBasis[k][1] == Dict([0, 1, 1, 0] => 2^0.5, [1, 1, 0, 0] => 2^0.5)
            @test transformedBasis[k][2] == Dict([0, 1, 1, 0] => 2^0.5, [1, 1, 0, 0] => -2^0.5)
            @test transformedBasis[k][3] == Dict([0, 0, 1, 1] => 1.0)
            @test transformedBasis[k][4] == Dict([1, 0, 0, 1] => 1.0)
        else
            @test transformedBasis[k] == basisStates[k]
        end
    end
end

@testset "Expand basis" begin
    basisStates = BasisStates(2)
    expandedBasis = expandBasis(basisStates, 1)
    @test issetequal(keys(expandedBasis), ((0, 0), (1, 1), (1, -1), (2, 0), (2, 2), (2, -2), (3, 1), (3, -1), (4, 0)))
    @test issetequal(expandedBasis[(0, 0)], (Dict([0, 0, 0, 0] => 1.0),))
    @test issetequal(expandedBasis[(1, 1)], (Dict([1, 0, 0, 0] => 1.0), Dict([0, 0, 1, 0] => 1.0)))
    @test issetequal(expandedBasis[(1, -1)], (Dict([0, 1, 0, 0] => 1.0), Dict([0, 0, 0, 1] => 1.0)))
    @test issetequal(expandedBasis[(2, 2)], (Dict([1, 0, 1, 0] => 1.0),))
    @test issetequal(expandedBasis[(2, 0)], (Dict([1, 1, 0, 0] => 1.0), Dict([0, 0, 1, 1] => 1.0), Dict([1, 0, 0, 1] => 1.0), Dict([0, 1, 1, 0] => 1.0)))
    @test issetequal(expandedBasis[(2, -2)], (Dict([0, 1, 0, 1] => 1.0),))
    @test issetequal(expandedBasis[(3, 1)], (Dict([1, 0, 1, 1] => 1.0), Dict([1, 1, 1, 0] => 1.0)))
    @test issetequal(expandedBasis[(3, -1)], (Dict([0, 1, 1, 1] => 1.0), Dict([1, 1, 0, 1] => 1.0)))
    @test issetequal(expandedBasis[(4, 0)], (Dict([1, 1, 1, 1] => 1.0),))
end
"""
