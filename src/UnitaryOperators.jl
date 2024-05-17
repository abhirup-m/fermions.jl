function KondoUnitaries(alpha, num_entangled, sectors)
    unitaryOperatorList = []
    IOMposition = 2 + 2 * num_entangled
    for sector in sectors
        IOMposition += 1
        if sector == 'h'
            for i in 1:num_entangled
                push!(unitaryOperatorList, [0.25 * alpha, "n+-", (1, IOMposition, 1 + 2 * i)])
                push!(unitaryOperatorList, [-0.25 * alpha, "n+-", (2, IOMposition, 1 + 2 * i)])
                push!(unitaryOperatorList, [-0.25 * alpha, "n+-", (1, IOMposition + 1, 2 + 2 * i)])
                push!(unitaryOperatorList, [0.25 * alpha, "n+-", (2, IOMposition + 1, 2 + 2 * i)])
                push!(unitaryOperatorList, [0.5 * alpha, "+-+-", (2, 1, IOMposition, 2 + 2 * i)])
                push!(unitaryOperatorList, [0.5 * alpha, "+-+-", (1, 2, IOMposition + 1, 1 + 2 * i)])
            end
        else
            for i in 1:num_entangled
                push!(unitaryOperatorList, [0.25 * alpha, "n+-", (1, 1 + 2 * i, IOMposition)])
                push!(unitaryOperatorList, [-0.25 * alpha, "n+-", (2, 1 + 2 * i, IOMposition)])
                push!(unitaryOperatorList, [-0.25 * alpha, "n+-", (1, 2 + 2 * i, IOMposition + 1)])
                push!(unitaryOperatorList, [0.25 * alpha, "n+-", (2, 2 + 2 * i, IOMposition + 1)])
                push!(unitaryOperatorList, [0.5 * alpha, "+-+-", (1, 2, 2 + 2 * i, IOMposition)])
                push!(unitaryOperatorList, [0.5 * alpha, "+-+-", (2, 1, 1 + 2 * i, IOMposition + 1)])
            end
        end
    end
        
    return unitaryOperatorList
end