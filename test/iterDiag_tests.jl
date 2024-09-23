#=@testset "Organise Operator" begin=#
#=    members = [1, 2, 3, 4]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members) == (operator, members, 1)=#
#=    end=#
#==#
#=    members = [5, 6, 3, 4]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members) == (operator, [3, 4, 5, 6], 1)=#
#=    end=#
#==#
#=    members = [1, 2, 5, 6]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members) == (operator, members, 1)=#
#=    end=#
#==#
#=    members = [5, 1, 3, 6]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members) == (operator, [1, 3, 5, 6], 1)=#
#=    end=#
#==#
#=    members = [1, 5, 2, 3]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-", "nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members) == (operator[[1, 3, 4, 2]], [1, 2, 3, 5], -1)=#
#=    end=#
#=    for bits in Iterators.product("+-nh", "+-nh", "nh", "+-")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members, newMembers) == (operator[[1, 3, 4, 2]], [1, 2, 3, 5], -1)=#
#=    end=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-", "+-")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members, newMembers) == (operator[[1, 3, 4, 2]], [1, 2, 3, 5], 1)=#
#=    end=#
#=    for bits in Iterators.product("+-nh", "+-nh", "nh", "nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members, newMembers) == (operator[[1, 3, 4, 2]], [1, 2, 3, 5], 1)=#
#=    end=#
#==#
#=    members = [1, 2, 5, 3]=#
#=    newMembers = [5, 6]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members, newMembers) == (operator[[1, 2, 4, 3]], [1, 2, 3, 5],=#
#=                                                                  ifelse(operator[4] ∈ "+-", -1, 1))=#
#=    end=#
#==#
#=    members = [1, 5, 2, 6]=#
#=    newMembers = [5, 6]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members, newMembers) == (operator[[1, 3, 2, 4]], [1, 2, 5, 6],=#
#=                                                                  ifelse(operator[3] ∈ "+-", -1, 1))=#
#=    end=#
#==#
#=    members = [1, 5, 6, 2]=#
#=    newMembers = [5, 6]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members, newMembers) == (operator[[1, 4, 2, 3]], [1, 2, 5, 6], 1)=#
#=    end=#
#==#
#=    members = [5, 1, 6, 2]=#
#=    newMembers = [5, 6]=#
#=    for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=        operator = join(bits)=#
#=        @test OrganiseOperator(operator, members, newMembers) == (operator[[2, 4, 1, 3]], [1, 2, 5, 6], =#
#=                                                                  ifelse(operator[2] ∈ "+-", -1, 1))=#
#=    end=#
#==#
#=    for oldPosition in 1:4=#
#=        members = [5,6,5,6]=#
#=        members[oldPosition] = 1=#
#=        newMembers = [5, 6]=#
#=        for bits in Iterators.product("+-nh", "+-nh", "+-nh", "+-nh")=#
#=            operator = join(bits)=#
#=            @test OrganiseOperator(operator, members, newMembers) == (operator, members, 1)=#
#=        end=#
#=    end=#
#=end=#


@testset "Iter Diag" begin
    totalSites = 3
    initSites = 1
    kondoJ = 1.
    hop_t = 0.9
    field = 1e-5 # to avoid degeneracies

    function getHamFlow(initSites::Int64, totalSites::Int64, hop_t::Float64, kondoJ::Float64)
        hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[]
        initHam = Tuple{String, Vector{Int64}, Float64}[]
        for site in 1:(initSites-1)
            push!(initHam, ("+-",  [1 + 2 * site, 3 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
            push!(initHam, ("+-",  [3 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
            push!(initHam, ("+-",  [2 + 2 * site, 4 + 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
            push!(initHam, ("+-",  [4 + 2 * site, 2 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}
        end
        for site in 1:initSites
            push!(initHam, ("nn",  [1, 1 + 2 * site], kondoJ/4)) # n_{d up, n_{0 up}
            push!(initHam, ("nn",  [1, 2 + 2 * site], -kondoJ/4)) # n_{d up, n_{0 down}
            push!(initHam, ("nn",  [2, 1 + 2 * site], -kondoJ/4)) # n_{d down, n_{0 up}
            push!(initHam, ("nn",  [2, 2 + 2 * site], kondoJ/4)) # n_{d down, n_{0 down}
            push!(initHam, ("+-+-",  [1, 2, 2 + 2 * site, 1 + 2 * site], kondoJ/2)) # S_d^+ S_0^-
            push!(initHam, ("+-+-",  [2, 1, 1 + 2 * site, 2 + 2 * site], kondoJ/2)) # S_d^- S_0^+
        end
        push!(initHam, ("n",  [1], field))
        push!(initHam, ("n",  [2], -field))

        push!(hamFlow, initHam)

        for site in initSites+1:totalSites
            newTerm = []
            push!(newTerm, ("+-",  [1 + 2 * site, -1 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
            push!(newTerm, ("+-",  [-1 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
            push!(newTerm, ("+-",  [2 + 2 * site, 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
            push!(newTerm, ("+-",  [2 * site, 2 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}
            push!(newTerm, ("nn",  [1, 1 + 2 * site], kondoJ/4)) # n_{d up, n_{0 up}
            push!(newTerm, ("nn",  [1, 2 + 2 * site], -kondoJ/4)) # n_{d up, n_{0 down}
            push!(newTerm, ("nn",  [2, 1 + 2 * site], -kondoJ/4)) # n_{d down, n_{0 up}
            push!(newTerm, ("nn",  [2, 2 + 2 * site], kondoJ/4)) # n_{d down, n_{0 down}
            push!(newTerm, ("+-+-",  [1, 2, 2 + 2 * site, 1 + 2 * site], kondoJ/2)) # S_d^+ S_0^-
            push!(newTerm, ("+-+-",  [2, 1, 1 + 2 * site, 2 + 2 * site], kondoJ/2)) # S_d^- S_0^+

            push!(hamFlow, newTerm)
        end
        return hamFlow
    end

    hamFlow = getHamFlow(initSites, totalSites, hop_t, kondoJ)
    exactHamiltonian = vcat(hamFlow...)
    basis = BasisStates(2 + 2 * totalSites)
    F = eigen(OperatorMatrix(basis, exactHamiltonian))

    @showprogress for members in collect(Iterators.product(fill(1:2+2*totalSites, 4)...))
        count = 1
        testingCorrelations = Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}()
        for bits in Iterators.product(fill("+-", 4)...)
            testingCorrelations[string(count)] = [(join(bits), collect(members), 1.)]
            count += 1
        end

        savePaths, _ = IterDiag(hamFlow, 1100, correlationDefDict = testingCorrelations, silent=true)
        results = deserialize.(savePaths)
        for (k, v) in testingCorrelations
            exactResult = F.vectors[:, 1]' * OperatorMatrix(basis, v) * F.vectors[:, 1]
            @test isapprox(exactResult, results[end-1]["results"][k][end], atol=1e-10)
        end

        testingCorrelations = Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}()
        for bits in Iterators.product(fill("nh", 4)...)
            testingCorrelations[string(count)] = [(join(bits), collect(members), 1.)]
            count += 1
        end
        savePaths, _ = IterDiag(hamFlow, 1100, correlationDefDict = testingCorrelations, silent=true)
        results = deserialize.(savePaths)
        for (k, v) in testingCorrelations
            exactResult = F.vectors[:, 1]' * OperatorMatrix(basis, v) * F.vectors[:, 1]
            @test isapprox(exactResult, results[end-1]["results"][k][end], atol=1e-10)
        end

        testingCorrelations = Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}()
        for bits in Iterators.product(fill("+n", 4)...)
            testingCorrelations[string(count)] = [(join(bits), collect(members), 1.)]
            count += 1
        end
        savePaths, _ = IterDiag(hamFlow, 1100, correlationDefDict = testingCorrelations, silent=true)
        results = deserialize.(savePaths)
        for (k, v) in testingCorrelations
            exactResult = F.vectors[:, 1]' * OperatorMatrix(basis, v) * F.vectors[:, 1]
            @test isapprox(exactResult, results[end-1]["results"][k][end], atol=1e-10)
        end
    end
end
