# compensation factor -⟨S_d^+ S_0^-⟩, where S_0 = \sum_{k1k2}S_{k1k2} for 1CK
# and S_0=\sum_{l1, l2} S^(l1, l2)} = \sum_l \sum_{k1 k2} S^(l1, l2)}_{k1 k2} for 2CK. The sum
# is dynamic, running over the k-states that are presently within the emergent
# window.
function CompensationOperator(states::Vector{Int64})
    compensationOperator = Dict{Tuple{String,Vector{Int64}},Float64}()
    for (i1, i2) in Iterators.product(states, states)
        mergewith!(+, compensationOperator, Dict(("+-+-", [1, 2, 2 * i1 + 2, 2 * i2 + 1]) => -1.0))
        mergewith!(+, compensationOperator, Dict(("n+-", [1, 2 * i1 + 1, 2 * i2 + 1]) => -0.5))
        mergewith!(+, compensationOperator, Dict(("n+-", [1, 2 * i1 + 2, 2 * i2 + 2]) => 0.5))
    end
    return [tuple(k..., v) for (k, v) in compensationOperator]
end


function ChargeFlip(leftStates::Vector{Int64}, rightStates::Vector{Int64})
    chargeFlipOperator = Tuple{String,Vector{Int64},Float64}[]
    for (i1, i2) in Iterators.product(leftStates, rightStates)
        push!(chargeFlipOperator, ("++--", [2 * i1 + 1, 2 * i1 + 2, 2 * i2 + 2, 2 * i2 + 1], 1.0 / (length(leftStates) * length(rightStates))))
    end
    return chargeFlipOperator
end


function SpinFlip(leftStates::Vector{Int64}, rightStates::Vector{Int64})
    spinFlip = Tuple{String,Vector{Int64},Float64}[]
    for (i1, i2) in Iterators.product(leftStates, rightStates)
        push!(spinFlip, ("+-+-", [2 * i1 + 1, 2 * i1 + 2, 2 * i2 + 2, 2 * i2 + 1], 1.0 / (length(leftStates) * length(rightStates))))
    end
    return spinFlip
end


function DegreeOfCompensation(leftStates2CK::Vector{Vector{Int64}}, states1CK::Vector{Vector{Int64}}, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}}, statesPerChannel::Int64)
    compensationOperator2CK = CompensationOperator.([leftStates2CK[1] for _ in leftStates2CK])
    compensationOperator1CK = CompensationOperator.([states1CK[1] for _ in states1CK])

    compensation1CK = fermions.FastCorrelation.(stateFlow1CK, compensationOperator1CK) ./ statesPerChannel
    compensation2CK = fermions.FastCorrelation.(stateFlow2CK, compensationOperator2CK) ./ statesPerChannel
    compensation2CKSI = fermions.FastCorrelation.(stateFlow2CKSI, compensationOperator2CK) ./ statesPerChannel
    return compensation2CK, compensation2CKSI, compensation1CK
end


function ChargeCorrelation(leftStates2CK::Vector{Vector{Int64}}, rightStates2CK::Vector{Vector{Int64}}, states1CK::Vector{Vector{Int64}}, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    chargeOperator2CK = [vcat(ChargeFlip(state1, state1), ChargeFlip(state1, state2)) for (state1, state2) in zip(leftStates2CK, rightStates2CK)]
    chargeOperator1CK = ChargeFlip.(states1CK, states1CK)
    chargeResults2CK = fermions.FastCorrelation.(stateFlow2CK, chargeOperator2CK)
    chargeResults2CKSI = fermions.FastCorrelation.(stateFlow2CKSI, chargeOperator2CK)
    chargeResults1CK = fermions.FastCorrelation.(stateFlow1CK, chargeOperator1CK)
    return chargeResults2CK, chargeResults2CKSI, chargeResults1CK
end


function SpinCorrelation(leftStates2CK::Vector{Vector{Int64}}, rightStates2CK::Vector{Vector{Int64}}, states1CK::Vector{Vector{Int64}}, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    spin2CK = [vcat(SpinFlip(state1, state1), SpinFlip(state1, state2)) for (state1, state2) in zip(leftStates2CK, rightStates2CK)]
    spin1CK = SpinFlip.(states1CK, states1CK)
    spinResults2CK = fermions.FastCorrelation.(stateFlow2CK, spin2CK)
    spinResults2CKSI = fermions.FastCorrelation.(stateFlow2CKSI, spin2CK)
    spinResults1CK = fermions.FastCorrelation.(stateFlow1CK, spin1CK)
    return spinResults2CK, spinResults2CKSI, spinResults1CK
end


function ImpVne(stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    vneImp2CKResults = [fermions.vnEntropy(state, [1, 2]; schmidtGap=true) for state in stateFlow2CK]
    vneImp2CKSIResults = [fermions.vnEntropy(state, [1, 2]; schmidtGap=true) for state in stateFlow2CKSI]
    vneImp1CKResults = [fermions.vnEntropy(state, [1, 2]; schmidtGap=true) for state in stateFlow1CK]
    vneImp2CK = [r[1] for r in vneImp2CKResults]
    vneImp2CKSI = [r[1] for r in vneImp2CKSIResults]
    vneImp1CK = [r[1] for r in vneImp1CKResults]
    schmidtGap2CK = [r[2] for r in vneImp2CKResults]
    schmidtGap2CKSI = [r[2] for r in vneImp2CKSIResults]
    schmidtGap1CK = [r[2] for r in vneImp1CKResults]
    return [vneImp2CK, vneImp2CKSI, vneImp1CK], [schmidtGap2CK, schmidtGap2CKSI, schmidtGap1CK]
end


function SchmidtGapBath(states2CKExpanded::Vector{Vector{Int64}}, states1CKExpanded::Vector{Vector{Int64}}, statesPerChannel::Int64, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    schmidtGap_kin_2CK = [maximum(fetch.([Threads.@spawn fermions.vnEntropy(state, states2CKExpanded[1][2*i-1:2*i]) for i in 1:2*statesPerChannel])) for state in stateFlow2CK]
    schmidtGap_kin_2CKSI = [maximum(fetch.([Threads.@spawn fermions.vnEntropy(state, states2CKExpanded[1][2*i-1:2*i]) for i in 1:2*statesPerChannel])) for state in stateFlow2CKSI]
    schmidtGap_kin_1CK = [maximum(fetch.([Threads.@spawn fermions.vnEntropy(state, states1CKExpanded[1][2*i-1:2*i]) for i in 1:statesPerChannel])) for state in stateFlow1CK]
    return schmidtGap_kin_2CK, schmidtGap_kin_2CKSI, schmidtGap_kin_1CK
end


function mutInfoInside(states2CKExpanded::Vector{Vector{Int64}}, states1CKExpanded::Vector{Vector{Int64}}, statesPerChannel::Int64, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    mutInfd_kin_2CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (states2CKExpanded[j][2*i-1:2*i], [1, 2])) for i in 1:2*statesPerChannel]))
                       for (j, state) in enumerate(stateFlow2CK)]
    mutInfd_kin_2CKSI = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (states2CKExpanded[j][2*i-1:2*i], [1, 2])) for i in 1:2*statesPerChannel]))
                         for (j, state) in enumerate(stateFlow2CKSI)]
    mutInfd_kin_1CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (states1CKExpanded[j][2*i-1:2*i], [1, 2])) for i in 1:statesPerChannel]))
                       for (j, state) in enumerate(stateFlow1CK)]
    return mutInfd_kin_2CK, mutInfd_kin_2CKSI, mutInfd_kin_1CK
end


function mutInfoOutside(states2CKExpanded::Vector{Vector{Int64}}, states1CKExpanded::Vector{Vector{Int64}}, statesPerChannel::Int64, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    mutInfd_kout_2CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, ([1, 2], states2CKExpanded[j][end-2*i+1:end-2*i+2])) for i in 1:2*statesPerChannel]))
                        for (j, state) in enumerate(stateFlow2CK)]
    mutInfd_kout_2CKSI = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, ([1, 2], states2CKExpanded[j][end-2*i+1:end-2*i+2])) for i in 1:2*statesPerChannel]))
                          for (j, state) in enumerate(stateFlow2CKSI)]
    mutInfd_kout_1CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, ([1, 2], states1CKExpanded[j][end-2*i+1:end-2*i+2])) for i in 1:statesPerChannel]))
                        for (j, state) in enumerate(stateFlow1CK)]
    return mutInfd_kout_2CK, mutInfd_kout_2CKSI, mutInfd_kout_1CK
end


function mutInfoInOut(states2CKExpanded::Vector{Vector{Int64}}, states1CKExpanded::Vector{Vector{Int64}}, statesPerChannel::Int64, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    mutInf_kout_kin_2CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (states2CKExpanded[j][2*i-1:2*i], states2CKExpanded[j][end-2*k+1:end-2*k+2])) for i in 1:2*statesPerChannel for k in 1:2*statesPerChannel]))
                           for (j, state) in enumerate(stateFlow2CK)]
    mutInf_kout_kin_2CKSI = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (states2CKExpanded[j][2*i-1:2*i], states2CKExpanded[j][end-2*k+1:end-2*k+2])) for i in 1:2*statesPerChannel for k in 1:2*statesPerChannel]))
                             for (j, state) in enumerate(stateFlow2CKSI)]
    mutInf_kout_kin_1CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (states1CKExpanded[j][2*i-1:2*i], states1CKExpanded[j][end-2*k+1:end-2*k+2])) for i in 1:statesPerChannel for k in 1:statesPerChannel]))
                           for (j, state) in enumerate(stateFlow1CK)]
    return mutInf_kout_kin_2CK, mutInf_kout_kin_2CKSI, mutInf_kout_kin_1CK
end


function mutInfoInLR(leftStates2CKExpanded::Vector{Vector{Int64}}, rightStates2CKExpanded::Vector{Vector{Int64}}, statesPerChannel::Int64, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}})
    mutInf_kin_LR_2CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (leftStates2CKExpanded[1][2*i-1:2*i], rightStates2CKExpanded[1][2*k-1:2*k])) for i in 1:statesPerChannel for k in 1:statesPerChannel]))
                         for (_, state) in enumerate(stateFlow2CK)]
    mutInf_kin_LR_2CKSI = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (leftStates2CKExpanded[1][2*i-1:2*i], rightStates2CKExpanded[1][2*k-1:2*k])) for i in 1:statesPerChannel for k in 1:statesPerChannel]))
                           for (_, state) in enumerate(stateFlow2CKSI)]
    return mutInf_kin_LR_2CK, mutInf_kin_LR_2CKSI
end


function mutInfoOutLR(leftStates2CKExpanded::Vector{Vector{Int64}}, rightStates2CKExpanded::Vector{Vector{Int64}}, statesPerChannel::Int64, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}})
    mutInf_kout_LR_2CK = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (leftStates2CKExpanded[j][end-2*i+1:end-2*i+2], rightStates2CKExpanded[j][2*k-1:2*k])) for i in 1:statesPerChannel for k in 1:statesPerChannel]))
                          for (j, state) in enumerate(stateFlow2CK)]
    mutInf_kout_LR_2CKSI = [maximum(fetch.([Threads.@spawn fermions.mutInfo(state, (leftStates2CKExpanded[j][end-2*i+1:end-2*i+2], rightStates2CKExpanded[j][2*k-1:2*k])) for i in 1:statesPerChannel for k in 1:statesPerChannel]))
                            for (j, state) in enumerate(stateFlow2CKSI)]
    return mutInf_kout_LR_2CK, mutInf_kout_LR_2CKSI
end


function tripInfoMax(states2CK_A::Vector{Vector{Int64}}, states2CK_B::Vector{Vector{Int64}}, states1CK_A::Vector{Vector{Int64}}, states1CK_B::Vector{Vector{Int64}}, statesPerChannel::Int64, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    tripInfo_2CK = [minimum(fetch.([Threads.@spawn fermions.tripartiteInfo(state, ([1, 2], states2CK_A[j][2*i-1:2*i], states2CK_B[j][2*k-1:2*k])) for i in 1:statesPerChannel for k in 1:statesPerChannel]))
                    for (j, state) in enumerate(stateFlow2CK)]
    tripInfo_2CKSI = [minimum(fetch.([Threads.@spawn fermions.tripartiteInfo(state, ([1, 2], states2CK_A[j][2*i-1:2*i], states2CK_B[j][2*k-1:2*k])) for i in 1:statesPerChannel for k in 1:statesPerChannel]))
                      for (j, state) in enumerate(stateFlow2CKSI)]
    tripInfo_1CK = [minimum(fetch.([Threads.@spawn fermions.tripartiteInfo(state, ([1, 2], states1CK_A[j][2*i-1:2*i], states1CK_B[j][2*k-1:2*k])) for i in 1:statesPerChannel for k in 1:statesPerChannel if i ≠ k]))
                    for (j, state) in enumerate(stateFlow1CK)]
    return tripInfo_2CK, tripInfo_2CKSI, tripInfo_1CK
end


function tripInfoBlock(states2CK_A::Vector{Vector{Int64}}, states2CK_B::Vector{Vector{Int64}}, states1CK_A::Vector{Vector{Int64}}, states1CK_B::Vector{Vector{Int64}}, stateFlow2CK::Vector{Dict{BitVector,Float64}}, stateFlow2CKSI::Vector{Dict{BitVector,Float64}}, stateFlow1CK::Vector{Dict{BitVector,Float64}})
    @assert states2CK_A ≠ states2CK_B
    @assert states1CK_A ≠ states1CK_B
    tripInfo_2CK = [fermions.tripartiteInfo(state, ([1, 2], states2CK_A[j], states2CK_B[j]))
                    for (j, state) in enumerate(stateFlow2CK)]
    tripInfo_2CKSI = [fermions.tripartiteInfo(state, ([1, 2], states2CK_A[j], states2CK_B[j]))
                      for (j, state) in enumerate(stateFlow2CKSI)]
    tripInfo_1CK = [fermions.tripartiteInfo(state, ([1, 2], states1CK_A[j], states1CK_B[j]))
                    for (j, state) in enumerate(stateFlow1CK)]
    return tripInfo_2CK, tripInfo_2CKSI, tripInfo_1CK
end
