function HubbardDimerMatrix(eps, U, hop_t)
    hubbardDimerMatrix = Dict()
    hubbardDimerMatrix[(0, 0)] = [0;;]
    hubbardDimerMatrix[(1, 1)] = [eps[3] hop_t[1]; hop_t[1] eps[1]]
    hubbardDimerMatrix[(1, -1)] = [eps[4] hop_t[2]; hop_t[2] eps[2]]
    hubbardDimerMatrix[(2, 2)] = [eps[1] + eps[3];;]
    hubbardDimerMatrix[(2, 0)] = (
                                        [U[2] + sum(eps[3:4]) hop_t[1] -hop_t[2] 0;
                                         hop_t[1] sum(eps[[1,4]]) 0 hop_t[2];
                                         -hop_t[2] 0 sum(eps[[2,3]]) -hop_t[1];
                                         0 hop_t[2] -hop_t[1] sum(eps[[1,2]]) + U[1]]
                                       )
    hubbardDimerMatrix[(2, -2)] = [eps[2] + eps[4];;]
    hubbardDimerMatrix[(3, 1)] = [U[2] + sum(eps[[1,3,4]]) -hop_t[2]; -hop_t[2] U[1] + sum(eps[1:3])]
    hubbardDimerMatrix[(3, -1)] = [U[2] + sum(eps[2:4]) -hop_t[1]; -hop_t[1] U[1] + sum(eps[[1,2,4]])]
    hubbardDimerMatrix[(4, 0)] = [sum(eps) + sum(U);;]
    return hubbardDimerMatrix
end


function HubbardDimerOplist(eps, U, hop_t)
    operatorList = Dict{Tuple{String, Vector{Int64}}, Float64}()
    operatorList[("n", [1])] = eps[1]
    operatorList[("n", [2])] = eps[2]
    operatorList[("n", [3])] = eps[3]
    operatorList[("n", [4])] = eps[4]
    operatorList[("nn", [1, 2])] = U[1]
    operatorList[("nn", [3, 4])] = U[2]
    operatorList[("+-", [1, 3])] = hop_t[1]
    operatorList[("+-", [3, 1])] = hop_t[1]
    operatorList[("+-", [2, 4])] = hop_t[2]
    operatorList[("+-", [4, 2])] = hop_t[2]
    return operatorList
end
