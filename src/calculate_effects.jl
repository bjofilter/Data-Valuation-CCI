using DataFrames

include("ida.jl")
include("functions.jl")


function add_all_effects(est_effects, multips, std_error)
    effs = DataFrame(Effect=Float64[], Probability=Float64[], StdError=Float64[])
    for k in keys(est_effects)
        if get(est_effects, k, 0) in effs[:,1]
            n = get(est_effects, k, 0)
            effs[effs.Effect .== n, :Probability] .= effs[effs.Effect .== n, :Probability] .+ get(multips, k, 0)/sum(values(multips))
        else
            push!(effs, (get(est_effects, k, 0), get(multips, k, 0)/sum(values(multips)), get(std_error, k, 0)))
        end
    end
    return effs
end


function calculate_CPDAGs!(coalitions, partitions, alpha)
    CPDAGs = [calculate_CPDAG(element, partitions, alpha) for element in coalitions.Members]
    coalitions[!, :CPDAG] = CPDAGs
    return coalitions
end


function calculate_CPDAG(coalition, partitions, alpha)
    println(coalition)
    data = partitions[coalition]
    data = combine_data(data)
    data = DataFrame(data, :auto)
    est_CPDAG = CausalInference.pcalg(data, alpha, gausscitest)
    # est_CPDAG = ges(data; penalty=1.0)[1]
    # est_CPDAG = exactscorebased(data; penalty=1.0)
    # if size(simplecycles(est_CPDAG), 1) > 0
    #     println("CYCLE DETECTED!!!!!")
    # end
    return est_CPDAG
end


function calculate_effect_all_pairs(CPDAG, CPDAG_all_data, coalition, partitions)

    println("coalition: ", coalition)

    data = partitions[coalition]
    data = combine_data(data)
    data = DataFrame(data, :auto)

    effects_all_pairs = DataFrame(x=Int64[], y=Int64[])
    foreach(x -> push!(effects_all_pairs, x), Iterators.product(1:nv(CPDAG), 1:nv(CPDAG)))
    effects_all_pairs[!, :Effects] = [DataFrame() for _ in 1:size(effects_all_pairs, 1)]
    effects_all_pairs = dropmissing!(ifelse.(effects_all_pairs.x .== effects_all_pairs.y, missing, effects_all_pairs))

    for i in 1:size(effects_all_pairs, 1)
        x, y = effects_all_pairs[i, :].x, effects_all_pairs[i, :].y
        if mod(i, 100) == 0
            println(i, " / ", size(effects_all_pairs, 1))
        end
        est_effects, multips, std_error = multiplicities(CPDAG, x, y, data)
        est_effects_orig_CPDAG = multiplicities(CPDAG_all_data, x, y, data)[1]
        if est_effects == est_effects_orig_CPDAG
            all_effects = add_all_effects(est_effects, multips, std_error)
            effects_all_pairs[i, :].Effects = all_effects
        else
            
        end
    end

    return effects_all_pairs
end