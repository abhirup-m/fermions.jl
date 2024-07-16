function HubbardDimerMatrix(eps, U, hop_t)
    hubbardDimerMatrix = Dict()
    hubbardDimerMatrix[(0, 0)] = [0;;]
    hubbardDimerMatrix[(1, 1)] = [eps[3] -hop_t[1]; -hop_t[1] eps[1]]
    hubbardDimerMatrix[(1, -1)] = [eps[4] -hop_t[2]; -hop_t[2] eps[2]]
    hubbardDimerMatrix[(2, 2)] = [eps[1] + eps[3];;]
    hubbardDimerMatrix[(2, 0)] = (
        [U[2]+sum(eps[3:4]) -hop_t[1] hop_t[2] 0;
        -hop_t[1] sum(eps[[1, 4]]) 0 -hop_t[2];
        hop_t[2] 0 sum(eps[[2, 3]]) hop_t[1];
        0 -hop_t[2] hop_t[1] sum(eps[[1, 2]])+U[1]]
    )
    hubbardDimerMatrix[(2, -2)] = [eps[2] + eps[4];;]
    hubbardDimerMatrix[(3, 1)] = [U[2]+sum(eps[[1, 3, 4]]) hop_t[2]; hop_t[2] U[1]+sum(eps[1:3])]
    hubbardDimerMatrix[(3, -1)] = [U[2]+sum(eps[2:4]) hop_t[1]; hop_t[1] U[1]+sum(eps[[1, 2, 4]])]
    hubbardDimerMatrix[(4, 0)] = [sum(eps) + sum(U);;]
    return hubbardDimerMatrix
end


function HubbardDimerSpecFunc(eps, U, hop_t, freqArr, broadening)
    Δ = (U^2 + 16 * hop_t^2)^0.5
    Egs = -U / 2 - Δ / 2

    a1 = 0.5 * √((Δ - U) / Δ)
    a2 = 2 * hop_t / √(Δ * (Δ - U))

    # basis = [|σ0>, |0σ>]
    exc_n_eq_1 = [a1, a2]

    eigPlus_n_1 = [1 / √2, 1 / √2]
    eigMinus_n_1 = [1 / √2, -1 / √2]
    En1 = [eps - hop_t, eps + hop_t]

    # basis = [|σ2>, |2σ>]
    exc_n_eq_3 = [a1, -a2]
    eigPlus_n_3 = [1 / √2, 1 / √2]
    eigMinus_n_3 = [1 / √2, -1 / √2]
    En3 = [eps + hop_t, eps - hop_t]

    A = broadening .* (
        sum(eigPlus_n_1 .* exc_n_eq_1)^2 ./ ((freqArr .- Egs .+ En1[1]) .^ 2 .+ broadening^2)
        .+
        sum(eigMinus_n_1 .* exc_n_eq_1)^2 ./ ((freqArr .- Egs .+ En1[2]) .^ 2 .+ broadening^2)
        .+
        sum(eigPlus_n_3 .* exc_n_eq_3)^2 ./ ((freqArr .+ Egs .- En3[1]) .^ 2 .+ broadening^2)
        .+
        sum(eigMinus_n_3 .* exc_n_eq_3)^2 ./ ((freqArr .+ Egs .- En3[2]) .^ 2 .+ broadening^2)
    )

    return A
end


function HubbardDimerOplist(eps, U, hop_t)
    operatorList = Tuple{String,Vector{Int64},Float64}[]
    if eps[1] != 0
        push!(operatorList, ("n", [1], eps[1]))
        push!(operatorList, ("n", [2], eps[2]))
        push!(operatorList, ("n", [3], eps[3]))
        push!(operatorList, ("n", [4], eps[4]))
    end
    if U[1] ≠ 0
        push!(operatorList, ("nn", [1, 2], U[1]))
        push!(operatorList, ("nn", [3, 4], U[2]))
    end
    if hop_t[1] ≠ 0
        push!(operatorList, ("+-", [1, 3], -hop_t[1]))
        push!(operatorList, ("+-", [3, 1], -hop_t[1]))
        push!(operatorList, ("+-", [2, 4], -hop_t[2]))
        push!(operatorList, ("+-", [4, 2], -hop_t[2]))
    end
    return operatorList
end


function dimerHamiltonian(U, t)
    hopping = Dict(
        ("+-", [1, 3]) => -t,
        ("+-", [2, 4]) => -t,
        ("+-", [3, 1]) => -t,
        ("+-", [4, 2]) => -t,
    )
    correlation = Dict(
        ("n", [1]) => -U / 2,
        ("n", [2]) => -U / 2,
        ("nn", [1, 2]) => U,
        ("n", [3]) => -U / 2,
        ("n", [4]) => -U / 2,
        ("nn", [3, 4]) => U,
    )
    return mergewith(+, hopping, correlation)
end


function dimerAdditionalHamiltonian(U, t, index)
    upPosition = 2 * index - 1
    hopping = Dict(
        ("+-", [upPosition, upPosition - 2]) => -t,
        ("+-", [upPosition + 1, upPosition - 1]) => -t,
        ("+-", [upPosition - 2, upPosition]) => -t,
        ("+-", [upPosition - 1, upPosition + 1]) => -t,
    )

    correlation = Dict(
        ("n", [upPosition]) => -U / 2,
        ("n", [upPosition + 1]) => -U / 2,
        ("nn", [upPosition, upPosition + 1]) => U
    )
    return mergewith(+, hopping, correlation)
end
