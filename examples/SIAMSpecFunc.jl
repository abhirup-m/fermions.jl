using fermions, Plots, Measures
include("../src/iterDiag.jl")

realSites = 200
maxSize = 2000
U = 8.
t = 1.

hamFlow = [
           [
           ("n", [1], -U/2),
           ("n", [2], -U/2),
           ("nn", [1, 2], U),
           ("+-",  [1, 3], -t),
           ("+-",  [3, 1], -t),
           ("+-",  [2, 4], -t),
           ("+-",  [4, 2], -t),
          ]
          ]
for site in 2:1:realSites
    newTerm = [
               ("+-",  [2 * site - 1, 2 * site + 1], -t),
               ("+-",  [2 * site + 1, 2 * site - 1], -t),
               ("+-",  [2 * site, 2 * site + 2], -t),
               ("+-",  [2 * site + 2, 2 * site], -t),
               #=("+-",  [2 * site + 1, 2 * site + 3], -t),=#
               #=("+-",  [2 * site + 3, 2 * site + 1], -t),=#
               #=("+-",  [2 * site + 2, 2 * site + 4], -t),=#
               #=("+-",  [2 * site + 4, 2 * site + 2], -t),=#
              ]
    push!(hamFlow, newTerm)
end

savePaths, _, specFuncOperators = IterDiag(hamFlow, maxSize;
                             symmetries=Char['N', 'S'],
                             # symmetries=Char['S'],
                             #=symmetries=Char['N'],=#
                             magzReq=(m, N) -> -2 ≤ m ≤ 3,
                             occReq=(x, N) -> div(N, 2) - 5 ≤ x ≤ div(N, 2) + 5,
                             #=corrMagzReq=(m, N) -> m == ifelse(isodd(div(N, 2)), 1, 0),=#
                             #=corrOccReq=(x, N) -> x == div(N, 2),=#
                             specFuncDefDict=Dict("create" => ("+", [1]), "destroy" => ("-", [1])),
                            )
#=specFuncOperators = Dict("destroy" => specFuncOperators["imp"][1], "create" => specFuncOperators["imp"][2])=#
freqValues = collect(-6:0.01:6)
@time specFunc = IterSpecFunc(savePaths, specFuncOperators, freqValues, 0.05;
                        occReq=(x,N) -> x == div(N,2),
                        excOccReq=(x,N) -> abs(x - div(N,2)) == 1,
                       )
p = plot(freqValues, specFunc)
savefig(p, "specfunc.pdf")
display(p)
