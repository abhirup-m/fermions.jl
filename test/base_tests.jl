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


@testset "RoundTo" begin
    @test RoundTo(1, 1e-10) == 1
    @test RoundTo(1 + 1e-11, 1e-10) == 1
    @test RoundTo(1e-11, 1e-10) == 0
    @test RoundTo(1e-10, 1e-10) == 1e-10
    @test RoundTo(-1e-10, 1e-10) == -1e-10
    @test RoundTo(0, 1e-10) == 0
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
    totalOperator = vcat(allOperators...)
    allStates = [Dict(k => coeffs2[i] * v for (k, v) in dict) for (i, dict) in enumerate(basis)]
    totalState = mergewith(+, allStates...)
    totalOutgoingState = ApplyOperator(totalOperator, totalState)
    pieceWiseAddedState = mergewith(+, [ApplyOperator(operator, state) for operator in allOperators for state in allStates]...)

    @test keys(totalOutgoingState) == keys(pieceWiseAddedState)
    for key in keys(totalOutgoingState)
        @test totalOutgoingState[key] ≈ pieceWiseAddedState[key]
    end

    # quantitative checking through basis states
    state = Dict(BitVector([0, 0]) => 0.1)
    @testset "state = [0, 0], operator = $op" for op in ["+", "-", "n", "h"]
        oplist = [(op, [1], 0.5)]
        if op == "+"
            @test ApplyOperator(oplist, state) == Dict(BitVector([1, 0]) => 0.05)
        elseif op == "h"
            @test ApplyOperator(oplist, state) == Dict(BitVector([0, 0]) => 0.05)
        else
            @test isempty(ApplyOperator(oplist, state))
        end
        oplist = [(op, [2], 0.5)]
        if op == "+"
            @test ApplyOperator(oplist, state) == Dict(BitVector([0, 1]) => 0.05)
        elseif op == "h"
            @test ApplyOperator(oplist, state) == Dict(BitVector([0, 0]) => 0.05)
        else
            @test isempty(ApplyOperator(oplist, state))
        end
        @testset for op2 in ["+", "-", "n", "h"]
            oplist = [(op * op2, [1, 2], 0.5)]
            if occursin("-", oplist[1][1]) || occursin("n", oplist[1][1])
                @test isempty(ApplyOperator(oplist, state))
            elseif oplist[1][1] == "++"
                @test ApplyOperator(oplist, state) == Dict(BitVector([1, 1]) => 0.05)
            elseif oplist[1][1] == "+h"
                @test ApplyOperator(oplist, state) == Dict(BitVector([1, 0]) => 0.05)
            elseif oplist[1][1] == "h+"
                @test ApplyOperator(oplist, state) == Dict(BitVector([0, 1]) => 0.05)
            elseif oplist[1][1] == "hh"
                @test ApplyOperator(oplist, state) == Dict(BitVector([0, 0]) => 0.05)
            end
        end
    end

    state = Dict(BitVector([1, 1]) => 0.1)
    @testset "state = [1, 1], operator = $op" for op in ["+", "-", "n", "h"]
        oplist = [(op, [1], 0.5)]
        if op == "-"
            @test ApplyOperator(oplist, state) == Dict(BitVector([0, 1]) => 0.05)
        elseif op == "n"
            @test ApplyOperator(oplist, state) == Dict(BitVector([1, 1]) => 0.05)
        else
            @test isempty(ApplyOperator(oplist, state))
        end
        oplist = [(op, [2], 0.5)]
        if op == "-"
            @test ApplyOperator(oplist, state) == Dict(BitVector([1, 0]) => -0.05)
        elseif op == "n"
            @test ApplyOperator(oplist, state) == Dict(BitVector([1, 1]) => 0.05)
        else
            @test isempty(ApplyOperator(oplist, state))
        end
        @testset "state = [1, 1], operator2 = $op2" for op2 in ["+", "-", "n", "h"]
            oplist = [(op * op2, [1, 2], 0.5)]
            if occursin("+", op * op2) || occursin("h", op * op2)
                @test isempty(ApplyOperator(oplist, state))
            elseif op * op2 == "--"
                @test ApplyOperator(oplist, state) == Dict(BitVector([0, 0]) => -0.05)
            elseif op * op2 == "-n"
                @test ApplyOperator(oplist, state) == Dict(BitVector([0, 1]) => 0.05)
            elseif op * op2 == "n-"
                @test ApplyOperator(oplist, state) == Dict(BitVector([1, 0]) => -0.05)
            elseif op * op2 == "nn"
                @test ApplyOperator(oplist, state) == Dict(BitVector([1, 1]) => 0.05)
            end
        end
    end
end

@testset "OperatorMatrix" begin
    eps = 1.281
    hop_t = 0.284
    U = 3.132
    operatorList = HubbardDimerOplist(eps, U, hop_t)
    @testset "sector=$((n, m))" for (n, m) in [(0, 0), (1, 1), (1, -1), (2, 2), (2, 0), (2, -2), (3, 1), (3, -1), (4, 0)]
        basisStates = BasisStates(4; totOccReq=n, magzReq=m)
        if (n,m) == (2,0)
            @test Set([collect(keys(b))[1] for b in basisStates]) == Set([[0, 0, 1, 1], [0, 1, 1, 0], [1, 0, 0, 1], [1, 1, 0, 0]])
            basisStates = [Dict(BitVector([1, 0, 0, 1]) => 1.0), Dict(BitVector([0, 1, 1, 0]) => 1.0), Dict(BitVector([1, 1, 0, 0]) => 1.0), Dict(BitVector([0, 0, 1, 1]) => 1.0)]
        end
        computedMatrix = OperatorMatrix(basisStates, operatorList)
        comparisonMatrix = HubbardDimerMatrix(eps, U, hop_t)[(n, m)]
        @test comparisonMatrix ≈ computedMatrix
    end
end


@testset "StateOverlap" begin
    basis = BasisStates(4)
    for (b1, b2) in Iterators.product(basis, basis)
        @test StateOverlap(b1, b2) == ifelse(b1 == b2, 1, 0)
    end
    coeffs1 = rand(length(basis))
    totalState1 = mergewith(+, [Dict(k => coeffs1[i] * v for (k, v) in dict) for (i, dict) in enumerate(basis)]...)
    coeffs2 = rand(length(basis))
    totalState2 = mergewith(+, [Dict(k => coeffs2[i] * v for (k, v) in dict) for (i, dict) in enumerate(basis)]...)
    @test StateOverlap(totalState1, totalState2) ≈ sum(coeffs1 .* coeffs2)
end
