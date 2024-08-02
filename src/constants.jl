const qbitMatrices = Dict{Char, Matrix{Float64}}()
qbitMatrices['n'] = diagm([0, 1.])
qbitMatrices['h'] = diagm([1., 0])
qbitMatrices['+'] = diagm(-1 => [1.])
qbitMatrices['-'] = diagm(1 => [1.])
const sigmax = qbitMatrices['+'] + qbitMatrices['-']
const sigmay = im .* qbitMatrices['+'] - im .* qbitMatrices['-']
const sigmaz = qbitMatrices['h'] - qbitMatrices['n']

export qbitMatrices, sigmax, sigmay, sigmaz
