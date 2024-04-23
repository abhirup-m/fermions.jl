include("fermionise.jl")
len = 18
basis = BasisStates(len, totOccupancy=[trunc(Int, len / 2)])
ham_op_list = [("nn", 1.0, [i, i + 1]) for i in 1:len-1]
H = generalOperatorMatrix(basis, oplist)
eigvals, eigstates = getSpectrum(H)
gstateCorrelation(basis, eigvals, eigstates, oplist)
