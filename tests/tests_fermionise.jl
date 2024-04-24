using Test
include("../source/fermionise.jl")

@testset "Basis States" begin
    b1 = BasisStates(1)
    @test Set(keys(b1)) == Set([(0, 0), (1, 1)])
    @test b1[(0, 0)][1] == [0]
    @test b1[(1, 1)][1] == [1]
    b2 = BasisStates(2)
    @test Set(keys(b2)) == Set([(0, 0), (1, 1), (1, -1), (2, 0)])
    @test b2[(0, 0)][1] == [0, 0]
    @test b2[(1, 1)][1] == [1, 0]
    @test b2[(1, -1)][1] == [0, 1]
    @test b2[(2, 0)][1] == [1, 1]
    b3 = BasisStates(3)
    @test Set(keys(b3)) == Set([(0, 0), (1, 1), (1, -1), (2, 0), (2, 2), (3, 1)])
    println(typeof(b2[(0, 0)]))
    @test Set(b3[(0, 0)]) == Set([BitArray([0, 0, 0])])
    @test Set(b3[(1, 1)]) == Set([BitArray([1, 0, 0]), BitArray([0, 0, 1])])
    @test Set(b3[(1, -1)]) == Set([BitArray([0, 1, 0])])
    @test Set(b3[(2, 0)]) == Set([BitArray([1, 1, 0]), BitArray([0, 1, 1])])
    @test Set(b3[(2, 2)]) == Set([BitArray([1, 0, 1])])
    @test Set(b3[(3, 1)]) == Set([BitArray([1, 1, 1])])
end


@testset "Transform Bit" begin
    @test TransformBit(true, 'n') == (1, 1)
    @test TransformBit(true, 'h') == (1, 0)
    @test TransformBit(true, '+') == (1, 0)
    @test TransformBit(true, '-') == (0, 1)
    @test TransformBit(false, 'n') == (0, 0)
    @test TransformBit(false, 'h') == (0, 1)
    @test TransformBit(false, '+') == (1, 1)
    @test TransformBit(false, '-') == (0, 0)
end


@testset "Apply Operator on State" begin
    # checking linearity
    b4 = BasisStates(4)
    for key in [(0, 0), (1, 1), (1, -1), (2, 0), (2, 2), (2, -2), (3, 1), (3, -1), (4, 0)]
        stateDict = Dict(k => v for (k, v) in zip(b4[key], rand(length(b4[key]))))
        partialStates = [Dict(k => v) for (k, v) in stateDict]
        ops = ["n", "h", "+", "-"]
        combinedOpList = Tuple{String,Float64,Vector{Int64}}[]
        totalApplication = []
        for op3 in Iterators.product(ops, ops, ops, ops)
            oplist = [(join(op3, ""), 1.0, [1, 2, 3, 4])]
            push!(combinedOpList, (join(op3, ""), 1.0, [1, 2, 3, 4]))
            partialApplications = [applyOperatorOnState(state, oplist) for state in partialStates]
            completeApplication = applyOperatorOnState(stateDict, oplist)
            @test completeApplication == merge(+, partialApplications...)
            totalApplication = [totalApplication; [completeApplication]]
        end
        @test merge(+, totalApplication...) == applyOperatorOnState(stateDict, combinedOpList)
    end

    # quantitative checking through basis states
    state = Dict(BitVector([0, 0]) => 0.1)
    @testset "state = [0, 0], operator = $op" for op in ["+", "-", "n", "h"]
        oplist = [(op, 0.5, [1])]
        if op == "+"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 0]) => 0.05)
        elseif op == "h"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 0]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        oplist = [(op, 0.5, [2])]
        if op == "+"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 1]) => 0.05)
        elseif op == "h"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 0]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        @testset for op2 in ["+", "-", "n", "h"]
            oplist = [(op * op2, 0.5, [1, 2])]
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
        oplist = [(op, 0.5, [1])]
        if op == "-"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([0, 1]) => 0.05)
        elseif op == "n"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 1]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        oplist = [(op, 0.5, [2])]
        if op == "-"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 0]) => -0.05)
        elseif op == "n"
            @test applyOperatorOnState(state, oplist) == Dict(BitVector([1, 1]) => 0.05)
        else
            @test isempty(applyOperatorOnState(state, oplist))
        end
        @testset for op2 in ["+", "-", "n", "h"]
            oplist = [(op * op2, 0.5, [1, 2])]
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


@testset "General operator matrix" begin
    b2 = BasisStates(4)
    oplist = [
        ("+-+-", 0.5, [1, 2, 4, 3]),
        ("+-+-", 0.5, [2, 1, 3, 4]),
    ]
    matrix = generalOperatorMatrix(b2, oplist)
    display(matrix[(2, 0)])
end
