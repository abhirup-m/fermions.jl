using Random

const op_cd = ('+', '-')
const op_nh = ('n', 'h')

function TestModel(totalSites::Int64)
    chemPot = 0.1
    hop_t_nn = rand(totalSites-1)
    hop_t_nnn = rand(totalSites-2)
    hubbard_U_nn = rand(totalSites-1)
    hubbard_U_3 = rand(totalSites-2)

    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for i in 1:totalSites-1
        push!(hamiltonian, ("+-", [i, i + 1], -hop_t_nn[i]))
        push!(hamiltonian, ("+-", [i + 1, i], -hop_t_nn[i]))
        push!(hamiltonian, ("nn", [i, i + 1], hubbard_U_nn[i]))
        push!(hamiltonian, ("n", [i], chemPot))
        if i < totalSites-1
            push!(hamiltonian, ("+-", [i, i + 2], -hop_t_nnn[i]))
            push!(hamiltonian, ("+-", [i + 2, i], -hop_t_nnn[i]))
            push!(hamiltonian, ("nnn", [i, i + 1, i + 2], hubbard_U_3[i]))
        end
    end
    push!(hamiltonian, ("n", [totalSites], chemPot))
    return hamiltonian
end

function GetCorrelationOperators(totalSites)
    qubitOperators = ["+", "-", "n", "h"]
    operators = Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}()

    # 1-particle correlations
    for _ in 1:5
        indices = rand(Int, 2) .|> abs .|> n -> rem(n, totalSites) .+ 1
        operators[randstring()] = [(join([shuffle(qubitOperators)[1] for _ in indices]), indices, 1.)]
    end

    # 2-particle correlations
    for _ in 1:5
        indices = rand(Int, 4) .|> abs .|> n -> rem(n, totalSites) .+ 1
        operators[randstring()] = [(join([shuffle(qubitOperators)[1] for _ in indices]), indices, 1.)]
    end
    return operators
end


function GetVNEParty(totalSites)
    parties = Dict{String, Vector{Int64}}()
    for size in 1:minimum((totalSites - 1, 5))
        for _ in 1:5
            party = shuffle(1:totalSites)[1:size]
            parties[randstring()] = sort(party)
        end
    end
    return parties
end


function GetMIParties(totalSites)
    parties = Dict{String, NTuple{2, Vector{Int64}}}()
    for size1 in 1:minimum((totalSites - 1, 4))
        for size2 in 1:minimum((4, totalSites - size1))
            party1 = shuffle(1:totalSites)[1:size1]
            party2 = shuffle([i for i in 1:totalSites if i ∉ party1])[1:size2]
            parties[randstring()] = (sort(party1), sort(party2))
        end
    end
    return parties
end


function GetSpecFuncChoices(totalSites)
    objects = Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}[]
    for size in 1:2
        for run in 1:4
            left = shuffle(1:totalSites)[1:size]
            right = shuffle(1:totalSites)[1:size]
            operator = join(shuffle(['+', '-', 'n', 'h'])[1:size])
            object = Dict("create" => [(operator, left, 1.)],
                                         "destroy" => Dagger([(operator, right, 1.)]),                                 
                                        )
            push!(objects, object)
        end
    end
    return objects
end

@testset "IterDiag Hubbard" begin
    @testset for totalSites in 3:4
        hamiltonian = TestModel(totalSites)
        basis = BasisStates(totalSites)
        exactEnergies, exactEigvecs = Spectrum(hamiltonian, basis)

        indexPartitions = collect(1:1:totalSites)
        hamFlow = MinceHamiltonian(hamiltonian, indexPartitions)
        operators = GetCorrelationOperators(totalSites)
        vneParties = GetVNEParty(totalSites)
        MIParties = GetMIParties(totalSites)
        savePaths, results, exitCode = IterDiag(hamFlow, 5000; correlationDefDict=copy(operators), vneDefDict=copy(vneParties), mutInfoDefDict=copy(MIParties), silent=true)
        energies = deserialize(savePaths[end-1])["eigVals"]

        @test length(energies) == length(exactEnergies)
        @testset for i in eachindex(exactEnergies)
            @test isapprox(exactEnergies[i], energies[i], atol=1e-11)
        end
        @testset for (name, operator) in operators
            exactCorr = GenCorrelation(exactEigvecs[1], operator)
            @test isapprox(exactCorr, results[name], atol=1e-11)
        end
        @testset for (name, party) in vneParties
            exactVNE = VonNEntropy(exactEigvecs[1], party)
            @test isapprox(exactVNE, results[name], atol=1e-11)
        end
        @testset for (name, parties) in MIParties
            exactMI = MutInfo(exactEigvecs[1], parties)
            @test isapprox(exactMI, results[name], atol=1e-11)
        end
    end
end

@testset "IterDiag Hubbard, symmetries" begin
    chemPot = 0.01
    @testset for totalSites in 3:4
        hamiltonian = TestModel(totalSites)

        basis = BasisStates(totalSites)
        exactEnergies, exactEigvecs = Spectrum(hamiltonian, basis)

        indexPartitions = collect(1:1:totalSites)
        hamFlow = MinceHamiltonian(hamiltonian, indexPartitions)
        operators = GetCorrelationOperators(totalSites)
        vneParties = GetVNEParty(totalSites)
        MIParties = GetMIParties(totalSites)
        savePaths, results, exitCode = IterDiag(hamFlow, 5000; correlationDefDict=copy(operators), vneDefDict=copy(vneParties), mutInfoDefDict=copy(MIParties), silent=true, symmetries=['N'])
        energies = deserialize(savePaths[end-1])["eigVals"]
        quantumNos = deserialize(savePaths[end-1])["quantumNos"]

        @test length(energies) == length(exactEnergies)
        totalNumOperator = [("n", [i], 1.) for i in 1:totalSites]
        @testset for i in eachindex(exactEnergies)
            @test isapprox(exactEnergies[i], energies[i], atol=1e-11)
            @test quantumNos[i][1] == round(GenCorrelation(exactEigvecs[i], totalNumOperator), digits=10)
        end
        @testset for (name, operator) in operators
            exactCorr = GenCorrelation(exactEigvecs[1], operator)
            @test isapprox(exactCorr, results[name], atol=1e-11)
        end
        @testset for (name, party) in vneParties
            exactVNE = VonNEntropy(exactEigvecs[1], party)
            @test isapprox(exactVNE, results[name], atol=1e-11)
        end
        @testset for (name, parties) in MIParties
            exactMI = MutInfo(exactEigvecs[1], parties)
            @test isapprox(exactMI, results[name], atol=1e-11)
        end
    end
end


@testset "Spectral Function" begin
    freqValues = collect(-10:0.01:10)
    standDev = 0.01
    @testset for totalSites in 6:9
        hamiltonian = TestModel(totalSites)
        specFuncDefDict = [Dict("create" => [("+", [totalSites], 1.)],
                                "destroy" => [("-", [totalSites], 1.)]
                               )
                          ]
        basis = BasisStates(totalSites)
        exactEnergies, exactEigvecs = Spectrum(hamiltonian, basis)

        indexPartitions = collect(1:1:totalSites)
        hamFlow = MinceHamiltonian(hamiltonian, indexPartitions)
        for specFuncOperator in specFuncDefDict
            savePaths, resultsDict, specFuncOperators = IterDiag(hamFlow, 5000; silent=true, specFuncDefDict=specFuncOperator)
            specFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, standDev)
            exactSpecFunc = SpecFunc(exactEnergies, exactEigvecs, specFuncOperator, freqValues, basis, standDev, normalise=false)
            error = (specFunc .- exactSpecFunc) .|> abs |> maximum
            @test isapprox(error, 0., atol=1e-10)
        end
    end
end
