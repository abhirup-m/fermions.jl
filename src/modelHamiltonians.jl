function KondoModel(
        numBathSites::Int64,
        hop_t::Float64,
        kondoJ::Float64;
        globalField::Float64=0.,
        couplingTolerance::Float64=1e-15,
    )
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]

    # intra-bath hopping
    if abs(hop_t) > couplingTolerance
        for site in 1:(numBathSites-1)
            push!(hamiltonian, ("+-",  [1 + 2 * site, 3 + 2 * site], -hop_t)) # c^†_{j,up} c_{j+1,up}
            push!(hamiltonian, ("+-",  [3 + 2 * site, 1 + 2 * site], -hop_t)) # c^†_{j+1,up} c_{j,up}
            push!(hamiltonian, ("+-",  [2 + 2 * site, 4 + 2 * site], -hop_t)) # c^†_{j,dn} c_{j+1,dn}
            push!(hamiltonian, ("+-",  [4 + 2 * site, 2 + 2 * site], -hop_t)) # c^†_{j+1,dn} c_{j,dn}
        end
    end

    # kondo terms
    if abs(kondoJ) > couplingTolerance
        push!(hamiltonian, ("nn",  [1, 3], kondoJ/4)) # n_{d up, n_{0 up}
        push!(hamiltonian, ("nn",  [1, 4], -kondoJ/4)) # n_{d up, n_{0 down}
        push!(hamiltonian, ("nn",  [2, 3], -kondoJ/4)) # n_{d down, n_{0 up}
        push!(hamiltonian, ("nn",  [2, 4], kondoJ/4)) # n_{d down, n_{0 down}
        push!(hamiltonian, ("+-+-",  [1, 2, 4, 3], kondoJ/2)) # S_d^+ S_0^-
        push!(hamiltonian, ("+-+-",  [2, 1, 3, 4], kondoJ/2)) # S_d^- S_0^+
    end

    # global magnetic field (to lift any trivial degeneracy)
    if abs(globalField) > couplingTolerance
        for site in 0:numBathSites
            push!(hamiltonian, ("n",  [1 + 2 * site], globalField))
            push!(hamiltonian, ("n",  [2 + 2 * site], -globalField))
        end
    end

    @assert !isempty(hamiltonian) "Hamiltonian is empty!"

    return hamiltonian
end
export KondoModel


function KondoModel(
        dispersion::Vector{Float64},
        kondoJ::Float64;
        globalField::Float64=0.,
        cavityIndices::Vector{Int64}=Int64[],
        couplingTolerance::Float64=1e-15,
    )
    numBathSites = length(dispersion)
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]

    # kinetic energy
    for site in 1:numBathSites
        if abs(dispersion[site]) < couplingTolerance
            continue
        end
        push!(hamiltonian, ("n",  [1 + 2 * site], dispersion[site])) # up spin
        push!(hamiltonian, ("n",  [2 + 2 * site], dispersion[site])) # down spin
    end

    # kondo terms
    if abs(kondoJ) > couplingTolerance
        for indices in Iterators.product(1:numBathSites, 1:numBathSites)
            if any(∈(cavityIndices), indices)
                continue
            end
            up1, up2 = 2 .* indices .+ 1
            down1, down2 = (up1, up2) .+ 1
            push!(hamiltonian, ("n+-",  [1, up1, up2], kondoJ_indices / 4)) # n_{d up, n_{0 up}
            push!(hamiltonian, ("n+-",  [1, down1, down2], -kondoJ_indices / 4)) # n_{d up, n_{0 down}
            push!(hamiltonian, ("n+-",  [2, up1, up2], -kondoJ_indices / 4)) # n_{d down, n_{0 up}
            push!(hamiltonian, ("n+-",  [2, down1, down2], kondoJ_indices / 4)) # n_{d down, n_{0 down}
            push!(hamiltonian, ("+-+-",  [1, 2, down1, up2], kondoJ_indices / 2)) # S_d^+ S_0^-
            push!(hamiltonian, ("+-+-",  [2, 1, up1, down2], kondoJ_indices / 2)) # S_d^- S_0^+
        end
    end

    # global magnetic field (to lift any trivial degeneracy)
    if abs(globalField) > couplingTolerance
        for site in 0:numBathSites
            push!(hamiltonian, ("n",  [1 + 2 * site], globalField))
            push!(hamiltonian, ("n",  [2 + 2 * site], -globalField))
        end
    end

    @assert !isempty(hamiltonian) "Hamiltonian is empty!"

    return hamiltonian
end
export KondoModel


function KondoModel(
        dispersion::Vector{Float64},
        kondoJ::Matrix{Float64};
        globalField::Float64=0.,
        couplingTolerance::Float64=1e-15,
        cavityIndices::Vector{Int64}=Int64[],
    )
    numBathSites = length(dispersion)
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]

    # kinetic energy
    for site in 1:numBathSites
        if abs(dispersion[site]) < couplingTolerance
            continue
        end
        push!(hamiltonian, ("n",  [1 + 2 * site], dispersion[site])) # up spin
        push!(hamiltonian, ("n",  [2 + 2 * site], dispersion[site])) # down spin
    end

    # kondo terms
    for indices in Iterators.product(1:numBathSites, 1:numBathSites)
        if any(∈(cavityIndices), indices)
            kondoJ_indices = 0
        else
            kondoJ_indices = kondoJ[indices...]
        end
        if abs(kondoJ_indices) < couplingTolerance
            continue
        end
        up1, up2 = 2 .* indices .+ 1
        down1, down2 = (up1, up2) .+ 1
        push!(hamiltonian, ("n+-",  [1, up1, up2], kondoJ_indices / 4)) # n_{d up, n_{0 up}
        push!(hamiltonian, ("n+-",  [1, down1, down2], -kondoJ_indices / 4)) # n_{d up, n_{0 down}
        push!(hamiltonian, ("n+-",  [2, up1, up2], -kondoJ_indices / 4)) # n_{d down, n_{0 up}
        push!(hamiltonian, ("n+-",  [2, down1, down2], kondoJ_indices / 4)) # n_{d down, n_{0 down}
        push!(hamiltonian, ("+-+-",  [1, 2, down1, up2], kondoJ_indices / 2)) # S_d^+ S_0^-
        push!(hamiltonian, ("+-+-",  [2, 1, up1, down2], kondoJ_indices / 2)) # S_d^- S_0^+
    end

    # global magnetic field (to lift any trivial degeneracy)
    if abs(globalField) > couplingTolerance
        for site in 0:numBathSites
            push!(hamiltonian, ("n",  [1 + 2 * site], globalField))
            push!(hamiltonian, ("n",  [2 + 2 * site], -globalField))
        end
    end

    @assert !isempty(hamiltonian) "Hamiltonian is empty!"

    return hamiltonian
end
export KondoModel


function KondoModel(
        dispersion::Vector{Float64},
        kondoJ::Float64,
        bathInt::Float64;
        intLegs::Int64=4,
        globalField::Float64=0.,
        couplingTolerance::Float64=1e-15,
        cavityIndices::Vector{Int64}=Int64[],
    )
    numBathSites = length(dispersion)
    hamiltonian = KondoModel(dispersion, kondoJ; globalField=globalField, 
                             cavityIndices=cavityIndices, couplingTolerance=couplingTolerance)

    if abs(bathInt) > couplingTolerance
        for indices in Iterators.product(repeat([1:numBathSites], 4)...)
            if length(unique(indices)) > intLegs
                continue
            end
            up1, up2, up3, up4 = 2 .* indices .+ 1
            down1, down2, down3, down4 = (up1, up2, up3, up4) .+ 1
            push!(hamiltonian, ("+-+-",  [up1, up2, up3, up4], -bathInt / 2)) # 
            push!(hamiltonian, ("+-+-",  [down1, down2, down3, down4], -bathInt / 2)) # 
            push!(hamiltonian, ("+-+-",  [up1, up2, down3, down4], bathInt / 2)) # 
            push!(hamiltonian, ("+-+-",  [down1, down2, up3, up4], bathInt / 2)) # 
        end
    end

    @assert !isempty(hamiltonian) "Hamiltonian is empty!"

    return hamiltonian
end
export KondoModel


function KondoModel(
        dispersion::Vector{Float64},
        kondoJ::Matrix{Float64},
        k_indices::Vector{Int64},
        bathIntFunc::Function;
        bathIntLegs::Int64=4,
        globalField::Float64=0.,
        couplingTolerance::Float64=1e-15,
        cavityIndices::Vector{Int64}=Int64[],
    )
    @assert length(k_indices) == length(dispersion)

    numBathSites = length(dispersion)
    hamiltonian = KondoModel(dispersion, kondoJ; globalField=globalField, 
                             cavityIndices=cavityIndices, couplingTolerance=couplingTolerance)
    for indices in Iterators.product(repeat([1:numBathSites], 4)...)
        bathIntVal = bathIntFunc([k_indices[i] for i in indices])
        if length(unique(indices)) > bathIntLegs || abs(bathIntVal) < couplingTolerance
            continue
        end
        up1, up2, up3, up4 = 2 .* indices .+ 1
        down1, down2, down3, down4 = (up1, up2, up3, up4) .+ 1
        push!(hamiltonian, ("+-+-",  [up1, up2, up3, up4], -bathIntVal / 2)) # 
        push!(hamiltonian, ("+-+-",  [down1, down2, down3, down4], -bathIntVal / 2)) # 
        push!(hamiltonian, ("+-+-",  [up1, up2, down3, down4], bathIntVal / 2)) # 
        push!(hamiltonian, ("+-+-",  [down1, down2, up3, up4], bathIntVal / 2)) # 
    end

    @assert !isempty(hamiltonian) "Hamiltonian is empty!"

    return hamiltonian
end
export KondoModel


function SiamKSpace(
        dispersion::Vector{Float64},
        hybridisation::Float64,
        impOnsite::Float64,
        impCorr::Float64;
        globalField::Float64=0.,
        cavityIndices::Vector{Int64}=Int64[],
    )
    numBathSites = length(dispersion)
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]

    # kinetic energy
    for site in 1:numBathSites
        push!(hamiltonian, ("n",  [1 + 2 * site], dispersion[site])) # up spin
        push!(hamiltonian, ("n",  [2 + 2 * site], dispersion[site])) # down spin
    end

    push!(hamiltonian, ("n",  [1], impOnsite)) # Ed nup
    push!(hamiltonian, ("n",  [2], impOnsite)) # Ed ndown
    push!(hamiltonian, ("nn",  [1, 2], impCorr)) # U nup ndown

    # hybridisation
    for site in 1:numBathSites
        if site ∈(cavityIndices)
            continue
        end
        if abs(hybridisation) < 1e-15
            continue
        end
        up = 2 * site + 1
        down = up + 1
        push!(hamiltonian, ("+-",  [1, up], hybridisation)) 
        push!(hamiltonian, ("+-",  [up, 1], hybridisation))
        push!(hamiltonian, ("+-",  [2, down], hybridisation))
        push!(hamiltonian, ("+-",  [down, 2], hybridisation))
    end

    # global magnetic field (to lift any trivial degeneracy)
    if globalField ≠ 0
        for site in 0:numBathSites
            push!(hamiltonian, ("n",  [1 + 2 * site], globalField))
            push!(hamiltonian, ("n",  [2 + 2 * site], -globalField))
        end
    end

    return hamiltonian
end
export SiamKSpace


function SiamKondoKSpace(
        dispersion::Vector{Float64},
        hybridisation::Float64,
        kondoJ::Float64,
        impOnsite::Float64,
        impCorr::Float64;
        globalField::Float64=0.,
        cavityIndices::Vector{Int64}=Int64[],
    )
    hamiltonian = SiamKSpace(dispersion, hybridisation, impOnsite, impCorr; globalField=globalField, cavityIndices=cavityIndices)

    numBathSites = length(dispersion)
    # kondo terms
    for indices in Iterators.product(1:numBathSites, 1:numBathSites)
        if any(∈(cavityIndices), indices)
            kondoJ_indices = 0
            continue
        else
            kondoJ_indices = kondoJ
        end
        if abs(kondoJ) < 1e-15
            continue
        end
        up1, up2 = 2 .* indices .+ 1
        down1, down2 = (up1, up2) .+ 1
        push!(hamiltonian, ("n+-",  [1, up1, up2], kondoJ_indices / 4)) # n_{d up, n_{0 up}
        push!(hamiltonian, ("n+-",  [1, down1, down2], -kondoJ_indices / 4)) # n_{d up, n_{0 down}
        push!(hamiltonian, ("n+-",  [2, up1, up2], -kondoJ_indices / 4)) # n_{d down, n_{0 up}
        push!(hamiltonian, ("n+-",  [2, down1, down2], kondoJ_indices / 4)) # n_{d down, n_{0 down}
        push!(hamiltonian, ("+-+-",  [1, 2, down1, up2], kondoJ_indices / 2)) # S_d^+ S_0^-
        push!(hamiltonian, ("+-+-",  [2, 1, up1, down2], kondoJ_indices / 2)) # S_d^- S_0^+
    end
    return hamiltonian

    
end
export SiamKondoKSpace


function Dispersion(
        numStates::Int64, 
        lowerCutOff::Float64,
        upperCutoff::Float64,
        discretisation::String;
        phSymmetry::Bool=true,
    )

    @assert abs(lowerCutOff) ≤ abs(upperCutoff)
    @assert discretisation == "log" || discretisation == "lin"
    if discretisation == "log"
        @assert lowerCutOff > 0
    end

    dispersion = zeros(numStates)
    if discretisation == "log"
        if phSymmetry
            if numStates % 2 == 0
                dispersion[1:2:end] .= 10. .^ range(log10(lowerCutOff), stop=log10(upperCutoff), length=div(numStates, 2))
                dispersion[2:2:end] .= -1 .* dispersion[1:2:end]
            else
                dispersion[3:2:end] .= 10. .^ range(log10(lowerCutOff), stop=log10(upperCutoff), length=div(numStates - 1, 2))
                dispersion[2:2:end] .= -1 .* dispersion[3:2:end]
            end
        else
            if numStates % 2 == 0
                dispersion .= 10. .^ range(log10(lowerCutOff), stop=log10(upperCutoff), length=numStates)
            else
                dispersion[2:end] .= 10. .^ range(log10(lowerCutOff), stop=log10(upperCutoff), length=numStates-1)
            end
        end
    else
        if phSymmetry
            if numStates % 2 == 0
                dispersion[1:2:end] = range(abs(lowerCutOff), stop=abs(upperCutoff), length=div(numStates, 2)) 
                dispersion[2:2:end] .= -1 .* dispersion[1:2:end]
            else
                @assert dispersion[1] == 0
                dispersion[1:2:end] = range(abs(lowerCutOff), stop=abs(upperCutoff), length=div(numStates+1, 2)) 
                dispersion[2:2:end] .= -1 .* dispersion[3:2:end]
            end
        else
            dispersion = collect(range(abs(lowerCutOff), stop=abs(upperCutoff), length=numStates))
        end
    end
    return dispersion
end
export Dispersion
