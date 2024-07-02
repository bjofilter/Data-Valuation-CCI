using Statistics
using Plots
using GLM
using Statistics
using StatsPlots

include("activelearning.jl")
include("causal_functions.jl")

function multiplicities(G, x, y, df)
    n = nv(G)
    memo = Dict{Set, BigInt}()
    fmemo = zeros(BigInt, n)
    effects = Dict{Vector{Int64}, Float64}()
    std_error = Dict{Vector{Int64}, Float64}()
    mults = Dict{Vector{Int64}, BigInt}()
    true_eff = Dict{Vector{Int64}, Float64}()
    x_y_connected = true


    # get directed parents of x
    bp = Vector{Int64}()
    for p in inneighbors(G, x)
        if !has_edge(G, x, p)
            push!(bp, p)
        end
    end
    
    U = copy(G)
    U.ne = 0
    for i = 1:n
        filter!(j->has_edge(G, j, i), U.fadjlist[i])
        filter!(j->has_edge(G, i, j), U.badjlist[i])
        U.ne += length(U.fadjlist[i])
    end
    base = 1
    xc = undef
    for component in connected_components(G)
        if x in component
            if !(y in component)
                x_y_connected = false
            end
        end
    end
    for component in connected_components(U)
        if !(x in component)
            cc = induced_subgraph(U, component)
            if !ischordal(cc[1])
                println("Undirected connected components are NOT chordal...Abort")
                println("Are you sure the graph is a CPDAG?")
                # is there anything more clever than just returning?
                return
            end
            base *= count(cc, memo, fmemo)
        else
            xc = component
        end
    end
    cc = induced_subgraph(U, xc)
    H = cc[1]
    mp = cc[2]
    ix = -1
    for i in 1:length(mp)
        if mp[i] == x
            ix = i
            break
        end
    end
    # est_cov = cov(Matrix(df))
    for parents in cliques(H, inneighbors(H, ix))
        resp = Term(Meta.parse("x" * string(y)))
        pr = ["x" * string(x)]
        pa = [x]
        for z in parents
            push!(pr, "x" * string(mp[z]))
            push!(pa, mp[z])
        end
        for z in bp
            push!(pr, "x" * string(z))
            push!(pa, z)
        end
        pr = map(x -> Meta.parse(x), pr)
        pred = Tuple(Term.(pr))
        mults[parents] = base * count_intsize(cc, ix, parents, memo, fmemo)
        if y in pa || !x_y_connected
            effects[parents] = 0
            std_error[parents] = 0
            true_eff[parents] = 0
        else
            # effects[parents] = reg_with_cov(est_cov, pa, y)[1]
            # model = lm(FormulaTerm(resp, pred), DataFrame(df, :auto))
            model = lm(FormulaTerm(resp, pred), df)
            effects[parents] = coef(model)[2]
            std = stderror(model)[2]
            if isnan(std)
                std_error[parents] = 0
            else
                std_error[parents] = std
            end
        end
    end
    return effects, mults, std_error
end


# function multiplicities(G, x, y, true_cov, df)
#         n = nv(G)
#         memo = Dict{Set, BigInt}()
#         fmemo = zeros(BigInt, n)
#         effects = Dict{Vector{Int64}, Float64}()
#         mults = Dict{Vector{Int64}, BigInt}()
#         true_eff = Dict{Vector{Int64}, Float64}()
    
    
#         # get directed parents of x
#         bp = Vector{Int64}()
#         for p in inneighbors(G, x)
#             if !has_edge(G, x, p)
#                 push!(bp, p)
#             end
#         end
        
#         U = copy(G)
#         U.ne = 0
#         for i = 1:n
#             filter!(j->has_edge(G, j, i), U.fadjlist[i])
#             filter!(j->has_edge(G, i, j), U.badjlist[i])
#             U.ne += length(U.fadjlist[i])
#         end
#         base = 1
#         xc = undef
#         for component in connected_components(U)
#             if !(x in component)
#                 cc = induced_subgraph(U, component)
#                 if !ischordal(cc[1])
#                     println("Undirected connected components are NOT chordal...Abort")
#                     println("Are you sure the graph is a CPDAG?")
#                     # is there anything more clever than just returning?
#                     return
#                 end
#                 base *= count(cc, memo, fmemo)
#             else
#                 xc = component
#             end
#         end
#         cc = induced_subgraph(U, xc)
#         H = cc[1]
#         mp = cc[2]
#         ix = -1
#         for i in 1:length(mp)
#             if mp[i] == x
#                 ix = i
#                 break
#             end
#         end
#         est_cov = cov(Matrix(df))
#         for parents in cliques(H, inneighbors(H, ix))
#             resp = Term(Meta.parse("x" * string(y)))
#             pr = ["x" * string(x)]
#             pa = [1]
#             for z in parents
#                 push!(pr, "x" * string(mp[z]))
#                 push!(pa, mp[z])
#             end
#             for z in bp
#                 push!(pr, "x" * string(z))
#                 push!(pa, z)
#             end
#             pr = map(x -> Meta.parse(x), pr)
#             pred = Tuple(Term.(pr))
#             mults[parents] = base * count_intsize(cc, ix, parents, memo, fmemo)
#             if y in pa
#                 effects[parents] = 0
#                 true_eff[parents] = 0
#             else
#                 effects[parents] = reg_with_cov(est_cov, pa, y)[1]
#                 true_eff[parents] = reg_with_cov(true_cov, pa, y)[1]
#             end
#         end
#     return effects, mults, true_eff
# end


function loc_ida(G, x, y, df)
    n = nv(G)
    effects = Dict{Vector{Int64}, Float64}()
    bp = Vector{Int64}()
    up = Vector{Int64}()
    for p in inneighbors(G, x)
        if !has_edge(G, x, p)
            push!(bp, p)
        else
            push!(up, p)
        end
    end
    mat = cov(Matrix(df))
    for parents in cliques(G, up)
        resp = Term(Meta.parse("x" * string(y)))
        pr = ["x" * string(x)]
        pa = [1]
        for z in parents
            push!(pr, "x" * string(z))
            push!(pa, z)
        end
        for z in bp
            push!(pr, "x" * string(z))
            push!(pa, z)
        end
        pr = map(x -> Meta.parse(x), pr)
        pred = Tuple(Term.(pr))
        if y in pa
            effects[parents] = 0
        else
            # effects[parents] = coef(lm(FormulaTerm(resp, pred), df))[2]
            effects[parents] = reg_with_cov(mat, pa, y)[1]
        end
    end
    return effects
end


function multiplicities_no_effects(G, x)
    n = nv(G)
    memo = Dict{Set, BigInt}()
    fmemo = zeros(BigInt, n)
    mults = Dict{Vector{Int64}, BigInt}()

    bp = Vector{Int64}()
    for p in inneighbors(G, x)
        if !has_edge(G, x, p)
            push!(bp, p)
        end
    end
    
    U = copy(G)
    U.ne = 0
    for i = 1:n
        filter!(j->has_edge(G, j, i), U.fadjlist[i])
        filter!(j->has_edge(G, i, j), U.badjlist[i])
        U.ne += length(U.fadjlist[i])
    end
    base = 1
    xc = undef
    for component in connected_components(U)
        if !(x in component)
            cc = induced_subgraph(U, component)
            # if !ischordal(cc[1])
            #     println("Undirected connected components are NOT chordal...Abort")
            #     println("Are you sure the graph is a CPDAG?")
            #     # is there anything more clever than just returning?
            #     return
            # end
            base *= count(cc, memo, fmemo)
        else
            xc = component
        end
    end
    cc = induced_subgraph(U, xc)
    H = cc[1]
    mp = cc[2]
    ix = -1
    for i in 1:length(mp)
        if mp[i] == x
            ix = i
            break
        end
    end
    for parents in cliques(H, inneighbors(H, ix))
        parentset = Vector{Int64}()
        for z in parents
            push!(parentset, mp[z])
        end
        for z in bp
            push!(parentset, z)
        end
        mults[parentset] = base * count_intsize(cc, ix, parents, memo, fmemo)
    end
    return mults
end



function loc_ida_mults(x, y, df, mults)
    effects = Dict{Vector{Int64}, Float64}()
    mat = cov(Matrix(df))

    for parents in mults.parentSet
        push!(parents, x)
        effects[parents] = reg_with_cov(mat, parents, y)[1]
    end
    return effects
end



function glob_ida_mults(x, y, df, mults)
    effects = Dict{Vector{Int64}, Float64}()
    multips = mults
    mat = cov(Matrix(df))

    for parents in mults.parentSet
        push!(parents, x)
        effects[parents] = reg_with_cov(mat, parents, y)[1]
    end
    return effects, multips
end