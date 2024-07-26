function HubbardDimerMatrix(eps, U, hop_t)
    ####### Calculations #######
    ###     H = -tΣ_σ(c^†_{1σ}c_{2σ} + h.c.) ϵΣ_i n_i + UΣ_i n_i↑ n_i↓
    ###     H|0,0⟩=0, H|2,2⟩=(4ϵ+2U)|2,2⟩, H|σ,σ⟩=2ϵ|σ,σ⟩,
    ###     H|σ,0⟩=ϵ|σ,0⟩-t|0,σ⟩, H|0,σ⟩=ϵ|0,σ⟩-t|σ,0⟩,
    ###     H|σ,2⟩=(3ϵ+U)|σ,2⟩+t|2,σ⟩, H|2,σ⟩=(3ϵ+U)|2,σ⟩+t|σ,2⟩,
    ###     H|2,0⟩=(2ϵ+U)|2,0⟩-t(|↑,↓⟩-|↓,↑⟩)
    ###     H|0,2⟩=(2ϵ+U)|0,2⟩-t(|↑,↓⟩-|↓,↑⟩)
    ###     H|σ,-σ⟩=2ϵ|σ,-σ⟩-tσ(|0,2⟩+|2,0⟩)
    hubbardDimerMatrix = Dict()

    # N = 0
    hubbardDimerMatrix[(0, 0)] = [0;;]

    # N= 4
    hubbardDimerMatrix[(4, 0)] = [4*eps + 2*U;;]
    
    # N = 1
    hubbardDimerMatrix[(1, 1)] = [eps -hop_t; -hop_t eps]
    hubbardDimerMatrix[(1, -1)] = [eps -hop_t; -hop_t eps]

    # N = 3
    hubbardDimerMatrix[(3, 1)] = [
                                  U+3*eps    hop_t;
                                  hop_t      U+3*eps
                                 ]
    hubbardDimerMatrix[(3, -1)] = [
                                   U+3*eps      hop_t;
                                   hop_t        U+3*eps
                                  ]

    # N = 2
    hubbardDimerMatrix[(2, 2)] = [2*eps;;]

    #   ↑↓      ↓↑      20      02
    hubbardDimerMatrix[(2, 0)] = (
        [
         2*eps       0           -hop_t        -hop_t;
         0           2*eps       hop_t         hop_t;
         -hop_t      hop_t       2*eps+U       0;
         -hop_t      hop_t       0             2*eps+U;
        ]
    )
    hubbardDimerMatrix[(2, -2)] = [2*eps;;]
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
    if eps != 0
        push!(operatorList, ("n", [1], eps))
        push!(operatorList, ("n", [2], eps))
        push!(operatorList, ("n", [3], eps))
        push!(operatorList, ("n", [4], eps))
    end
    if U ≠ 0
        push!(operatorList, ("nn", [1, 2], U))
        push!(operatorList, ("nn", [3, 4], U))
    end
    if hop_t ≠ 0
        push!(operatorList, ("+-", [1, 3], -hop_t))
        push!(operatorList, ("+-", [3, 1], -hop_t))
        push!(operatorList, ("+-", [2, 4], -hop_t))
        push!(operatorList, ("+-", [4, 2], -hop_t))
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
