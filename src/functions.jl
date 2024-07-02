using Combinatorics
using Distributions
using Graphs
using InformationGeometry
using Random


function assignweights(G)
    n = nv(G)
    wg = [Dict{Integer, AbstractFloat}() for i = 1:n]
    ud = Uniform(1, 2)
    for a in vertices(G), b in outneighbors(G, a)
        wg[a][b] = rand(ud)
    end
    return wg
end


function calculate_marginal_contributions(num_partitions, coalitions, minV)
    marginal_contributions = DataFrame(Member=Int64[], Coalition=Vector{Int64}[], Contribution=Float64[])
    combinations = [i for i in powerset(1:num_partitions, 1)]
    for i in 1:num_partitions
        V_i = coalitions[coalitions.Members .== [[i]], :].V[1]
        push!(marginal_contributions, (i, [], V_i - minV))
        for j in combinations
            if 5 > size(j)[1] && size(j)[1] > 0 && !(i in j)
                k = copy(j)
                sort!(append!(k, i))
                V_j = coalitions[coalitions.Members .== [k], :].V[1]
                push!(marginal_contributions, (i, j, V_j - V_i))
            end
        end
    end
    return marginal_contributions
end


function calculate_modified_shapley_values(x, max_sv, rho, minV, coalitions)
    i = x.Member
    m1 = coalitions[coalitions.Members .== [[i]], :].V[1] - minV
    m2 = -minV * (x.SV / max_sv)^rho
    modified_shapley_value = maximum([m1, m2]) + minV
    return modified_shapley_value
end


function calculate_reverse_KLS(distributions_all_data, distributions_partition)
    true_distributions = distributions_all_data[:, :Models]
    est_distributions = distributions_partition[:, :Models]
    dist = 0
    count = 0
    for i in eachindex(true_distributions)
        m1 = true_distributions[i]
        m2 = est_distributions[i]
        if m2 != 0
            kb = KullbackLeibler(m2, m1)
            if !isnan(kb)
             dist -= kb
             count = count + 1
            end
        end
    end

    return dist / size(true_distributions, 1)
end


function calculate_shapley_values(num_partitions, marginal_contributions)
    shapley_values = DataFrame(Member=Int64[], SV=Float64[])
    combinations = [i for i in powerset(1:num_partitions, 1)]
    for i in 1:num_partitions
        sv = 0
        for j in combinations
            if !(i in j)
                m = marginal_contributions[findall((marginal_contributions.Member .== [i]) .& (marginal_contributions.Coalition .== [j])), :].Contribution
                sv = sv + factorial(size(j)[1]) * factorial(num_partitions - size(j)[1] - 1) * m[1]
            end
        end
        sv = sv / factorial(num_partitions)
        push!(shapley_values, (i, sv))
    end
    return shapley_values
end



function combine_data(data)
    d = data[1]
    for i in 2:size(data, 1)
        d = vcat(d, data[i])
    end
    return d
end


function compute_distributions(effects)
    models = [get_distribution(row.Effects) for row in eachrow(effects)]
    effects[!, :Models] = models
    return effects
end


function computeSHDs(a, b)
    n_row = size(a, 1)
    dist = 0
    for i in 1:n_row
        a_bits = BitArray(a[i, :])
        b_bits = BitArray(b[i, :])
        for i in 1:length(a_bits)
            dist += a_bits[i] ⊻ b_bits[i]
        end
    end
    return dist/n_row^2
end


function computeSIDs(cp_a, cp_b, num_nodes)
    different_effects = 0
    for x in 1:num_nodes
        for y in 1:num_nodes
            if x != y && !contain_same_effects(cp_a, cp_b, x, y)
                different_effects += 1
            end
        end
    end
    return different_effects / (num_nodes^2 - num_nodes)
end


function contain_same_effects(cp_a, cp_b, x, y)

    multips_a = multiplicities_no_effects(cp_a, x)
    multips_b = multiplicities_no_effects(cp_b, x)
    parentSets_a = get_true_diff_parentSets(multips_a, cp_a, x, y)
    parentSets_b = get_true_diff_parentSets(multips_b, cp_b, x, y)

    for pset_a in parentSets_a.parentSet
        if valid_adj_set(pset_a, parentSets_b, cp_b, x, y)
            return true
        end
    end
    for pset_b in parentSets_b.parentSet
        if valid_adj_set(pset_b, parentSets_a, cp_a, x, y)
            return true
        end
    end

    return false
end


function get_coalitions(num_partitions)
    coalitions = DataFrame()
    combinations = [i for i in powerset(1:num_partitions, 1)]
    coalitions[!, :Members] = combinations
    return coalitions
end


function get_distribution(effects)
    if size(effects)[1] == 0
        model = 0
    else
        normal_dists = [Normal(row.Effect, row.StdError) for row in eachrow(effects)]
        priors = effects[:, :Probability]
        model = MixtureModel(normal_dists, priors)
    end
    return model
end


function makeDAG!(G, ts)
    n = nv(G)
    for a in 1:n, b in inneighbors(G, a)
        if ts[b] < ts[a]
            rem_edge!(G, a, b)
        end
    end
end


function randgraph(n, d, method)
    G = SimpleDiGraph(n)
    if method == "ba"
        D = Graphs.barabasi_albert!(Graphs.random_regular_digraph(Int64(round(n*d)), Int64(round(n*d*d))), n, Int64(round((n-1)*d)))
        for a = 1:n, b in neighbors(D, a)
            add_edge!(G, a, b)
            add_edge!(G, b, a)
        end
    elseif method == "kb"
        ud = Uniform(0, 1)
        for a = 1:n, b in 1:n
            if a < b && rand(ud) < d
                add_edge!(G, a, b)
                add_edge!(G, b, a)
            end
        end
    else
        return "no valid method"
    end
    ts = randperm(num_nodes)
    makeDAG!(G, ts)
    wg = assignweights(G)
    return G, wg, ts
end


function sampleData(wg, s, ts)
    n = size(wg, 1)
    nd = Normal()
    dt = rand(nd,s, n)
    its = zeros(Int64, n)
    for i in 1:n
        its[ts[i]] = i
    end
    for a in its, b in wg[a]
        dt[1:s, b.first] += dt[1:s, a] * b.second
    end
    return dt
end


function valid_adj_set(pset_a, parentSets_b, cp_b, x, y)

    for pset_b in parentSets_b.parentSet
        if is_adj_set(pset_a, pset_b, cp_b, x, y)
            return true
        end
    end
    return false
end