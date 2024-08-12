@testset "Iter Diag" begin
    realSites = 2
    maxSize = 2
    Ed = -1.
    t = 1.
    hamFlow = [
               [
               ("n", [1], Ed),
              ]
              ]
    for site in 1:realSites
        newTerm = [
                   ("+-",  [site, site + 1], -t),
                   ("+-",  [site + 1, site], -t),
                  ]
        push!(hamFlow, newTerm)
    end

    numSitesFlow = collect(1:1+realSites)
    savePaths = IterDiag(hamFlow, maxSize);
    files = [deserialize(path) for path in savePaths]
    results = [Dict("rotation" => f["rotation"], "eigVals" => f["eigVals"]) for f in files]
    @test results[1]["eigVals"] ≈ [Ed, 0]
    @test results[1]["rotation"] ≈ [[0 1]; [1 0]]
    F = eigen([[Ed 0 0 -t]; [0 Ed 0 0]; [0 0 0 0]; [-t 0 0 0]])
    @test results[2]["eigVals"] ≈ F.values[1:2]
    @test results[2]["rotation"] ≈ F.vectors[:, 1:2]
    F = eigen([[results[2]["eigVals"][1] 0 0 0];
             [0 results[2]["eigVals"][1] -t * F.vectors[:, 1][1] 0];
               [0 -t * F.vectors[:, 1][1] results[2]["eigVals"][2] 0];
               [0 0 0 results[2]["eigVals"][2]]
              ]
             )
    @test results[3]["eigVals"] ≈ F.values[1:2]
    @test results[3]["rotation"] ≈ F.vectors[:, 1:2]
end
