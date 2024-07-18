function kondoOperators(momIndices::NTuple{2,Integer}, multiplier::Float64)
    return [
        ("n+-", [1, momIndices[1], momIndices[2]], 0.25 * multiplier),
        ("n+-", [2, momIndices[1], momIndices[2]], -0.25 * multiplier),
        ("n+-", [1, 1 + momIndices[1], momIndices[2] + 1], -0.25 * multiplier),
        ("n+-", [2, 1 + momIndices[1], momIndices[2] + 1], 0.25 * multiplier),
        ("+-+-", [1, 2, 1 + momIndices[1], momIndices[2]], 0.5 * multiplier),
        ("+-+-", [2, 1, momIndices[1], momIndices[2] + 1], 0.5 * multiplier),
    ]
end


function unitaries2CK(alpha::Float64, num_entangled::Integer, sectors::String)
    # num_entangled indicates the number of lattice sites (including 
    # impurity site) already present within model before the IOM that
    # is about to be inserted.

    # labelling scheme (L(R)i = conduction electrons of left(right) channel)
    # d↑  d↓    L1↑   L1↓   R1↑   R1↓   L-IOM-P↑   L-IOM-P↓   L-IOM-H↑   L-IOM-H↓   R-IOM-P↑   R-IOM-P↓   R-IOM-H↑   R-IOM-H↓
    # 1   2     3     4     5     6     7          8          9          10         11         12         13         14
    # num_entangled = 3 in this example.

    @assert sectors ∈ ("p", "h", "ph")
    @assert (num_entangled - 1) % 2 == 0

    # all odd sites are the left channel sites.
    # the total number of channel sites is num_entangled-1,
    # -1 being the impurity site
    leftChannel = 1:2:(num_entangled-1)

    # all even sites are the right channel sites
    rightChannel = 2:2:(num_entangled-1)

    # first left channel IOM is inserted at num_entangled + 1, 
    # num_entangled being the total number of sites prior to insertion
    # and hence also the index of the last site before insertion.
    leftIOM = num_entangled + 1

    # first right channel IOM site is inserted just after the first
    # left channel IOM site.
    rightIOM = leftIOM + 1
    unitaryTerms = []
    if 'p' in sectors
        leftPTerms = vcat([kondoOperators((2 * i + 1, 2 * leftIOM + 1), alpha) for i in leftChannel]...)
        rightPTerms = vcat([kondoOperators((2 * i + 1, 2 * rightIOM + 1), alpha) for i in rightChannel]...)
        unitaryTerms = [unitaryTerms; leftPTerms; rightPTerms]

        # next left IOM will be two sites away from the previous leftIOM,
        # because there's also a right IOM site in between
        leftIOM += 2

        # next right IOM site is also two sites away, because the second
        # right IOM has now been inserted between them.
        rightIOM += 2
    end
    if 'h' in sectors
        leftHTerms = vcat([kondoOperators((2 * leftIOM + 1, 2 * i + 1), alpha) for i in leftChannel]...)
        rightHTerms = vcat([kondoOperators((2 * rightIOM + 1, 2 * i + 1), alpha) for i in rightChannel]...)
        unitaryTerms = [unitaryTerms; leftHTerms; rightHTerms]
    end
    return unitaryTerms
end


function unitaries1CK(alpha::Float64, num_entangled::Integer, sectors::String)
    # num_entangled indicates the number of lattice sites already present within
    # model before the IOM that is about to be inserted.

    # labelling scheme
    # d↑    d↓      k1↑     k1↓     k2↑     k2↓     k3↑     k3↓
    # 1     2       3       4       5       6       7       8 
    # num_entangled = 4 in this example.

    @assert sectors ∈ ("p", "h", "ph")
    IOMposition = 2 * num_entangled + 1
    innerStates = 2:num_entangled
    unitaryTerms = []
    if 'p' in sectors
        particleTerms = vcat([kondoOperators((2 * i - 1, IOMposition), alpha) for i in innerStates]...)
        unitaryTerms = [unitaryTerms; particleTerms]
        IOMposition += 2
    end
    if 'h' in sectors
        holeTerms = vcat([kondoOperators((IOMposition, 2 * i - 1), alpha) for i in innerStates]...)
        unitaryTerms = [unitaryTerms; holeTerms]
    end
    return unitaryTerms
end
