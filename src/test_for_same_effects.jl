using GraphPlot
using DataStructures



function get_effects_adj_sets(df, x, y, G)

    parents = inneighbors(G, x)
    cpDAG = cpdag(G)

    est_effects = DataFrame(Parents = String[], Effect=Float64[], Multip = Float64[], IDA=String[])
    true_effect = 0
        
    effects, multips = multiplicities(cpDAG, x, y, df)
    effects_loc = loc_ida(cpDAG, x, y, df)

    for k in keys(effects)
        if rep == 1
            push!(est_effects,(string(k), get(effects, k, 0), multips[k], "global IDA"))
            push!(est_effects,(string(k), get(effects, k, 0), 1, "local IDA"))
        else
            est_effects[est_effects.Parents .== string(k), :Effect] .+= get(effects, k, 0)
        end

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

    return est_effects, true_effect

end


function get_true_diff_parentSets(multips, cp, x, y)

    parentSets = collect(keys(multips))
    multips_new = copy(multips)
    pSets  = DataFrame(parentSet = Vector{Int64}[], equalSets = String[], multip = Float64[])

    for i in 1:length(parentSets)
        if parentSets[i] ∉ keys(multips_new)
            continue
        else
            for j in (i+1):length(parentSets)
                if parentSets[j] ∉ keys(multips_new)
                    continue
                elseif is_adj_set(parentSets[i], parentSets[j], cp, x, y)
                    if parentSets[i] in pSets.parentSet
                        pSets[pSets.parentSet .== [parentSets[i]], :equalSets] .= string(pSets[pSets.parentSet .== [parentSets[i]], :].equalSets[1], ", ", parentSets[j])
                        pSets[pSets.parentSet .== [parentSets[i]], :multip] .= pSets[pSets.parentSet .== [parentSets[i]], :multip][1] + multips[parentSets[j]]
                    else
                        push!(pSets,(parentSets[i], string(parentSets[j]), multips[parentSets[i]] + multips[parentSets[j]]))
                    end
                    multips_new[parentSets[i]] += multips[parentSets[j]]
                    delete!(multips_new, parentSets[j])
                end
            end
            if !(parentSets[i] in pSets.parentSet)
                push!(pSets,(parentSets[i], "none", multips[parentSets[i]]))
                #?
                delete!(multips_new, parentSets[i])
            end
        end
    end

    return pSets
end


function is_adj_set(a, b, cp, x, y)

    if y in a && y in b
        return true
    end

    B = adjustGraph(cp, x, b)

    if check_first(B, a, x, y)
        if check_second(B, a, x, y)
            return true
        else
            return false
        end
    else
       return false
    end
end


function adjustGraph(cp, x, parents)
    G = copy(cp)
    for i in all_neighbors(G, x)
        if i in parents
            rem_edge!(G, x, i)
        else
            rem_edge!(G, i, x)
        end
    end

    meek!(G)

    return G
end


function check_first(G, Z, x, y)
    A = copy(G)

    for i in outneighbors(A, x)
        if has_path(A, i, y)
            for z in Z
                if has_path(A, i, z)
                    return false
                end
            end
        end
    end
    return true
end


function check_second(G, Z, x, y)

    function visit(f, j)
        if !marked[[f, j]] && j != x
            enqueue!(q, [f, j])
            marked[[f, j]] = true
        end
    end

    n = nv(G)
    marked = Dict()
    q = Queue{Vector{Int64}}()

    g = copy(G)
    for v in outneighbors(G, x)
        if !(v in inneighbors(G, x)) && has_path(G, v, y)
            rem_edge!(g, x, v)
        end
    end

    for node in 1:n
        for neighbor in inneighbors(g, node)
            marked[[node, neighbor]] = false
        end
        for neighbor in outneighbors(g, node)
            marked[[node, neighbor]] = false
        end
    end
    

    for v in outneighbors(g, x)
        visit(x, v)
    end
    for v in inneighbors(g, x)
        visit(x, v)
    end

    while length(q) > 0
        f, k = dequeue!(q)

        if k == y
            return false
        end


        # f - k
        if !(k in Z) && f in inneighbors(g, k) && f in outneighbors(g, k)
            for node in outneighbors(g, k)
                visit(k, node)
            end
        end

        # f -> k
        if k in Z && f in inneighbors(g, k) && !(f in outneighbors(g, k))
            for node in inneighbors(g, k)
                visit(k, node)
            end
        end

        if !(k in Z) && f in inneighbors(g, k) && !(f in outneighbors(g, k))
            for node in outneighbors(g, k)
                if !(node in inneighbors(g, k))
                    visit(k, node)
                end
            end
        end

        # f <- k
        if !(k in Z) && !(f in inneighbors(g, k)) && f in outneighbors(g, k)
            for node in inneighbors(g, k)
                visit(k, node)
            end
            for node in outneighbors(g, k)
                visit(k, node)
            end
        end
    end

    return true
end


function effects_multip_for_parents(pSets, effects)
    est_effects = DataFrame(Parents = Vector{Int64}[], Effect=Float64[], Multip = Float64[], IDA=String[])

    for parents in eachrow(pSets)
        println(get(effects, parents.parentSet, 0))
        push!(est_effects,(parents.parentSet, get(effects, parents.parentSet, 0), parents.multip_glob, "global IDA"))
        push!(est_effects,(parents.parentSet, get(effects, parents.parentSet, 0), parents.multip_loc, "local IDA"))
    end

    return est_effects
end