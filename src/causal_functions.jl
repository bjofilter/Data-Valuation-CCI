using LinearSolve

# Given a DAG and two node x and y, compute the effect of x on y by multiplying the weights along each path x to y and then summing over
# all possible paths
function get_effect_from_dag(dag, x, y, wg, df, path_weight = 1, visited = [], effect = 0)
    if !has_path(dag, x, y)
        return 0
    else
        push!(visited, x)

        for a in outneighbors(dag, x)
            if has_path(dag, a, y)

                if a in keys(wg[x])
                    weight = wg[x][a]
                else
                   weight = wg[a][x] * var(df[!, string("x", a)]) / var(df[!, string("x", x)])
                end

                if a in visited
                    error("circle detected")
                elseif a == y
                    effect = effect + path_weight * weight
                else
                    effect = effect + get_effect_from_dag(dag, a, y, wg, df, path_weight * weight, copy(visited))
                end

            end
        end
    return effect
    end
end


function get_true_effects_from_true_dag(G, x, y, wg, parent_sets)
    effects = Dict()
    for parents in keys(parent_sets)
        effect = effect_from_true_dag(G, x, y, wg, parents)

        if effect in keys(effects)
            effects[effect] = effects[effect] + parent_sets[parents]
        else
            merge!(effects, Dict(effect => parent_sets[parents]))
        end
    end

    return effects
end


function effect_from_true_dag(dag, x, y, wg, parents, path_weight = 1, visited = [], effect = 0)
    if !connected(dag, x, y)
        return 0
    else
        push!(visited, x)

        for a in vcat(neighbors(dag, x), inneighbors(dag, x))
            if a in parents
                continue
            elseif connected(dag, a, y)

                if a in keys(wg[x])
                    weight = wg[x][a]
                else
                    weight = 1/wg[a][x]
                end

                if a in visited
                    continue
                elseif a == y
                    effect = effect + path_weight * weight
                else
                    effect = effect + get_effect_from_dag(dag, a, y, wg, path_weight * weight, copy(visited))
                end

            end
        end
    return effect
    end
end


# Get the effects of x on y for all DAGs (all possible extensions of one CPDAG)
function get_true_effects(dags, x, y, wg)
    true_effects = Dict()
    for i in 1:length(dags)
        true_effect = get_effect_from_dag(dags[i], x, y, wg)
        if true_effect in keys(true_effects)
            true_effects[true_effect] = true_effects[true_effect] + 1
        else
            merge!(true_effects, Dict(true_effect => 1))
        end
    end
    return true_effects
end


# Computes the effect of x on y via sampled data
function compute_only_true_reg(wg, ts, s, reps, x, y, alpha, parent_sets, dags)

    rep = 1
    true_effects = DataFrame(Effect=Float64[], Parents=String[])
    while rep <= reps
        dt = sampledata(wg, s, ts)
        df = DataFrame(dt, :auto)

        for parents in keys(parent_sets)
            resp = Term(Meta.parse("x" * string(y)))
            pr = ["x" * string(x)]
            for z in parents
                push!(pr, "x" * string(z))
            end
            pr = map(x -> Meta.parse(x), pr)
            pred = Tuple(Term.(pr))
            effect = coef(lm(FormulaTerm(resp, pred), df))[2]
            # push!(true_effects, (effect, string(parents)))
            for i in 1:parent_sets[parents]
                push!(true_effects, (effect, "effect reg"))
            end
        end

        for dag in dags
            effect_dag = get_effect_from_dag(dag, x, y, wg, df)
            push!(true_effects, (effect_dag, "effect dag"))
        end

        rep += 1
    end

    return true_effects

end


function compare_global_vs_local(wg, ts, s, reps, x, y, G)

    parents = inneighbors(G, x)
    cp = cpdag(G)

    rep = 1
    est_effects = DataFrame(Effect=Float64[], multip = Float64[], IDA=String[])
    true_effect = 0
    while rep <= reps
        dt = sampledata(wg, s, ts)
        df = DataFrame(dt, :auto)
        
        effects, multips = multiplicities(cp, x, y, df)
        effects_loc = loc_ida(cp, x, y, df)

        println(multips)

        for k in keys(effects)
            # push!(est_effects,(get(effects, k, 0), multips[k], "global IDA"))
            for i in 1:multips[k]
                push!(est_effects,(get(effects, k, 0), 1, "global IDA"))
            end
        end

        for k in keys(effects_loc)
            push!(est_effects,(get(effects_loc, k, 0), 1, "local IDA"))
        end

        resp = Term(Meta.parse("x" * string(y)))
        pr = ["x" * string(x)]
        for z in parents
            push!(pr, "x" * string(z))
        end
        pr = map(x -> Meta.parse(x), pr)
        pred = Tuple(Term.(pr))
        true_reg = coef(lm(FormulaTerm(resp, pred), df))[2]
        true_effect += true_reg

        rep += 1
    end

    true_effect /= reps

    return est_effects, true_effect

end


function connected(dag, a, b)
    for i in 1:nv(dag)
        if has_path(dag, i, a) && has_path(dag, i, b)
            return true
        end
    end
    return false
end


function percentage_in_range(est_effects, true_effect)
    num_global = sum(est_effects[est_effects.IDA .== "global IDA", :multip])
    num_local = sum(est_effects[est_effects.IDA .== "local IDA", :multip])

    range = abs(true_effect) * 0.01
    if range < 0.01
        range = 0.01
    end

    effects_in_range = est_effects[true_effect - range .< est_effects.Effect .< true_effect + range, :]

    num_global_in_range = sum(effects_in_range[effects_in_range.IDA .== "global IDA", :multip])
    num_local_in_range = sum(effects_in_range[effects_in_range.IDA .== "local IDA", :multip])

    percent_global = num_global_in_range / num_global
    percent_local = num_local_in_range / num_local
    
    return percent_global, percent_local
end


function precentage_exact(est_parent_sets, amount_true)
    percent_local = 1 / size(est_parent_sets)[1]

    percent_global = amount_true / sum(est_parent_sets.multip)

    return percent_global, percent_local
end


function get_multips(x, G)

    parents = inneighbors(G, x)
    cp = cpdag(G)

    parent_sets = DataFrame(Parents=Vector{Int64}[], multip = Float64[])
    amount_true = 0

    multips = multiplicities_no_effects(cp, x)

    for k in keys(multips)
        push!(parent_sets,(k, get(multips, k, 0)))
        if k == parents
            amount_true = get(multips, k, 0)
        end
    end

    return parent_sets, parents, amount_true
end

function get_effect_parents(pSets, G, x, y, wg, s, ts)

    effects = Dict{Vector{Int64}, Float64}()
    n = nv(G)

    dt = sampledata(wg, s, ts)
    df = DataFrame(dt, :auto)

    bp = Vector{Int64}()
    for p in inneighbors(G, x)
        if !has_edge(G, x, p)
            push!(bp, p)
        end
    end

    for parents in eachrow(pSets)
        resp = Term(Meta.parse("x" * string(y)))
        pr = ["x" * string(x)]
        for z in parents.parentSet
            push!(pr, "x" * string(z))
        end
        for z in bp
            push!(pr, "x" * string(z))
        end
        pr = map(x -> Meta.parse(x), pr)
        pred = Tuple(Term.(pr))
        effects[parents.parentSet] = coef(lm(FormulaTerm(resp, pred), df))[2]
    end
    return effects
end

function get_cov(G, ts, wg)
    n = nv(G)
    var_mat = zeros(n, n)
    its = zeros(Int64, n)
    for i in 1:n
        its[ts[i]] = i
    end
    for i in 1:n
        node = its[i]
        neigh = inneighbors(G, node)
        if (size(neigh, 1) == 0)
            var_mat[node, node] = 1
        else
            var_i = 1
            for neig_a in neigh
                for neig_b in neigh
                    var_i += wg[neig_a][node] * wg[neig_b][node] * var_mat[neig_a, neig_b]
                end
            end
            var_mat[node, node] = var_i
        end
        for j in 1:(i-1)
            node_b = its[j]
            var_ij = 0
            for neig_a in neigh
                var_ij += wg[neig_a][node] * var_mat[neig_a, node_b]
            end
            var_mat[node, node_b] = var_ij
            var_mat[node_b, node] = var_ij

        end




    end


    return var_mat
end

function get_cov_true_est(G, ts, wg)
    n = nv(G)
    nd = Normal()
    var_mat = zeros(n, n)
    var_mat_est = zeros(n, n)
    its = zeros(Int64, n)
    for i in 1:n
        its[ts[i]] = i
    end
    for i in 1:n
        node = its[i]
        neigh = inneighbors(G, node)
        if (size(neigh, 1) == 0)
            var_mat[node, node] = 1
            var_mat_est[node, node] = 1
        else
            var_i = 1
            var_i_est = 1
            for neig_a in neigh
                for neig_b in neigh
                    var_i += wg[neig_a][node] * wg[neig_b][node] * var_mat[neig_a, neig_b]
                    var_i_est += wg[neig_a][node] * wg[neig_b][node] * var_mat_est[neig_a, neig_b]
                end
            end
            err = rand(nd) * Float64(2.47)^-(n-i)
            var_mat[node, node] = var_i
            var_mat_est[node, node] = var_i_est + err
        end
        for j in 1:(i-1)
            node_b = its[j]
            neigh_b = inneighbors(G, node_b)
            var_ij = 0
            var_ij_est = 0
            for neig_a in neigh
                var_ij += wg[neig_a][node] * var_mat[neig_a, node_b]
                var_ij_est += wg[neig_a][node] * var_mat_est[neig_a, node_b]
            end
            var_mat[node, node_b] = var_ij
            var_mat[node_b, node] = var_ij
            var_mat_est[node, node_b] = var_ij_est
            var_mat_est[node_b, node] = var_ij_est

        end
    end
    return var_mat, var_mat_est
end

function reg_with_cov(mat, parents, y)
    parents = sort(parents)
    n = size(parents,1)
    A = zeros(n, n)
    for i in 1:n
        A[i, :] = mat[parents[i], parents]
    end
    b = mat[y, parents]

    prob = LinearProblem(A, b)
    coeff = solve(prob)
    return coeff
end

function get_std(effects, total=false)
    effects_for_dense = DataFrame(Effect=Float64[], IDA=String[])
    maxi = maximum(effects.Multip)
    while maxi > 10000
        effects.Multip = effects.Multip ./ sqrt(maxi)
        maxi = maximum(effects.Multip)
    end
    for j in 1:nrow(effects)
        for k in 1:round(effects[j, :Multip])
            push!(effects_for_dense,(effects[j, :Effect], "totel IDA"))
        end
    end
    std_e = std(effects_for_dense.Effect)

    return std_e
end