using GraphPlot
using CausalInference

include("calculate_effects.jl")
include("functions.jl")
include("test_for_same_effects.jl")

# Define variables
method_graph_construction = "kb"
alpha = 0.001
density = 0.9
num_nodes = 9
num_partitions = 5
s = 40000
rho = 0.9

# Create data
G, wg, ts = randgraph(num_nodes, density, method_graph_construction)
data = sampleData(wg, s, ts)

# Partition data into batches
batchsize = Int64(s / num_partitions)
partitions = [data[i:i+batchsize-1, :] for i in 1:batchsize:s]

# Combine batches to get all coalitions
coalitions = get_coalitions(num_partitions)

# Compute CPDAGS and adjacency matrixes
println("calculating CPDAGs...")
coalitions = calculate_CPDAGs!(coalitions, partitions, alpha)
CPDAG_all_data = last(coalitions.CPDAG)
adjacency_CPDAG_partitions = Matrix.(adjacency_matrix.(coalitions.CPDAG))
adjacency_CPDAG_all_data = last(adjacency_CPDAG_partitions)

# Calculate SHDs
SHDs = [computeSHDs(adjacency_CPDAG_all_data, adjacency_CPDAG_partitions[i]) for i in 1:size(coalitions, 1)]
coalitions[!, :SHD] = SHDs

# Calculate new SIDs
SIDs = [computeSIDs(coalitions.CPDAG[i], CPDAG_all_data, num_nodes) for i in 1:size(coalitions, 1)]
coalitions[!, :SID] = SIDs

# Calculate causal effects for all possible x-y pairs
println("calculate effects for all x-y pairs ...")
effect_distributions_partitions = [calculate_effect_all_pairs(row.CPDAG, CPDAG_all_data, row.Members, partitions) for row in eachrow(coalitions)]
# coalitions[!, :Effects] = effect_distributions_partitions
effect_distributions_all_data = last(effect_distributions_partitions)

# Compute distributions
println("computing distributions ...")
distributions_partitions = [compute_distributions(element) for element in effect_distributions_partitions]
# coalitions[!, :Distributions] = distributions_partitions
distributions_all_data = last(distributions_partitions)

# Calculate reverse KL-distances, averaged over all effects
reverse_KLs = [calculate_reverse_KLS(distributions_all_data, element) for element in distributions_partitions]
 coalitions[!, :NegReverseKL] = reverse_KLs

# Calculate coalition values
 coalitions[!, :V] = -coalitions[:, :SHD] +  coalitions[:, :NegReverseKL]
minV = minimum(coalitions[!, :V])
marginal_contributions = calculate_marginal_contributions(num_partitions, coalitions, minV)

# Calculate shapley values
shapley_values = calculate_shapley_values(num_partitions, marginal_contributions)
max_sv = maximum(shapley_values.SV)
modified_shapley_values = [calculate_modified_shapley_values(row, max_sv, rho, minV, coalitions) for row in eachrow(shapley_values)]
shapley_values[!, :MSV] = modified_shapley_values

# nodelabel = 1:num_nodes
# gplot(CPDAG_all_data, nodelabel=nodelabel)
# simplecycles(CPDAG_all_data)

gplot(dg, nodelabel=nodelabel)