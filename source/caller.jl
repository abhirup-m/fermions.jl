using BenchmarkTools
include("fermionise.jl")
len = 16
@btime _ = BasisStates(len, totOccupancy=[trunc(Int, len / 2)]);
basis = BasisStates(len, totOccupancy=[trunc(Int, len / 2)])
oplist = [("nn", 1.0, [i, i + 1]) for i in 1:len-1]
@btime _ = generalOperatorMatrix(basis, oplist);
H = generalOperatorMatrix(basis, oplist)
@btime _ = getSpectrum(H)
eigvals, eigstates = getSpectrum(H)
@btime gstateCorrelation(basis, eigvals, eigstates, oplist)
