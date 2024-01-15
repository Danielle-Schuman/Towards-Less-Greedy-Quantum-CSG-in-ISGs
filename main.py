import numpy as np
import random
import pickle
import time
import datetime

import utils
import plotting
from algorithm import Algorithm
from GCSQ import GCSQ
from jonas import Jonas
from danielle import Danielle
from n_split_GCSQ import n_split_GCSQ
from ours_iterative_exactly import ours_iterative_exactly
from ours_iterative_at_most import ours_iterative_at_most
from k_split_GCSQ_exactly import k_split_GCSQ_exactly
from k_split_GCSQ_at_most import k_split_GCSQ_at_most
from r_qubo import r_qubo
from r_qubo_iterative import r_qubo_iterative

def create_data_synthetic_test(graph_sizes, num_graphs_per_size, seed, mean=0.5):
    data = {}
    for n in graph_sizes:
        graphs = []
        for graph_num in range(num_graphs_per_size):
            graph = utils.generate_problem(n, mean=mean)
            graphs.append(graph)
        data[n] = graphs
    pickle.dump(data, open(f"data/tests/data_{graph_sizes}_{num_graphs_per_size}_{seed}.pkl", 'wb'))
    return data, True

# Setting the seed
seed = random.randint(0, 2**32 - 1)
random.seed(seed)
np.random.seed(seed)
print(f"Seed: {seed}")

# loading E.ON data
data, synthetic = pickle.load(open('data/data_new_20samples_4_28.pkl', 'rb')), False
# alternative: load synthetic data
# data, synthetic = create_data_synthetic_test([10, 15, 20], 25, seed)

graph_sizes = list(data.keys())  # is [4,6,8,10,12,14,16,18,20,22,24,26,28] for E.ON data
num_graphs_per_size = len(data[graph_sizes[0]])  # is 20 for E.ON data
num_graph_sizes = len(graph_sizes)

# TODO: Find sensible k's to try, possibly everything between 2 and n
k = 4
utils.qaoa_init(seed=seed)
algorithm_list = [GCSQ(seed=seed, num_graph_sizes=num_graph_sizes), Jonas(seed=seed, num_graph_sizes=num_graph_sizes), Danielle(seed=seed, num_graph_sizes=num_graph_sizes), n_split_GCSQ(seed=seed, num_graph_sizes=num_graph_sizes),
                  ours_iterative_exactly(seed=seed, num_graph_sizes=num_graph_sizes, k=k), ours_iterative_at_most(seed=seed, num_graph_sizes=num_graph_sizes, k=k), k_split_GCSQ_exactly(seed=seed, num_graph_sizes=num_graph_sizes, k=k), k_split_GCSQ_at_most(seed=seed, num_graph_sizes=num_graph_sizes, k=k),
                  r_qubo(seed=seed, num_graph_sizes=num_graph_sizes), r_qubo_iterative(seed=seed, num_graph_sizes=num_graph_sizes, k=k),
                  GCSQ(seed=seed, num_graph_sizes=num_graph_sizes, solver="qaoa"), Jonas(seed=seed, num_graph_sizes=num_graph_sizes, solver="qaoa"), Danielle(seed=seed, num_graph_sizes=num_graph_sizes, solver="qaoa"), n_split_GCSQ(seed=seed, num_graph_sizes=num_graph_sizes, solver="qaoa"),
                  ours_iterative_exactly(seed=seed, num_graph_sizes=num_graph_sizes, k=k, solver="qaoa"), ours_iterative_at_most(seed=seed, num_graph_sizes=num_graph_sizes, k=k, solver="qaoa"), k_split_GCSQ_exactly(seed=seed, num_graph_sizes=num_graph_sizes, k=k, solver="qaoa"), k_split_GCSQ_at_most(seed=seed, num_graph_sizes=num_graph_sizes, k=k, solver="qaoa"),
                  r_qubo(seed=seed, num_graph_sizes=num_graph_sizes, solver="qaoa"), r_qubo_iterative(seed=seed, num_graph_sizes=num_graph_sizes, k=k, solver="qaoa")]
results_dict = {}
run_id = str(datetime.datetime.now().date()) + '_' + str(datetime.datetime.now().time()).replace(':', '-')

for i, num_agents in enumerate(graph_sizes):
    print(f"\n\n\nTest for graphsize {num_agents}")
    results_dict[num_agents] = {}
    for graph_num in range(num_graphs_per_size):
        print(f"\n\n     Graph {graph_num}")
        graph = data[num_agents][graph_num]
        if synthetic:
            edges = graph
        else:
            edges = utils.transform(graph)
        for algorithm in algorithm_list:
            print(f"\n          Running {algorithm.name}...")
            start_time = time.time()
            coalitions = algorithm.solve(num_agents, edges)
            end_time = time.time()
            value = np.sum([utils.value(c, edges) for c in coalitions])
            total_time = end_time - start_time
            algorithm.data.append((coalitions, value, total_time))
            pickle.dump(algorithm.data, open(f"results/data_{algorithm.name}_{run_id}.pkl", 'wb'))
            results_dict[num_agents].setdefault(algorithm.name, []).append(value)
            print(f"          Coalition structure value for {algorithm.name}: {value}   -   Time: {total_time}")
            algorithm.values_sums[i] += value
    print(f"Results for Test with graph size {num_agents}: {results_dict[num_agents]}\n")
print("Done running tests, plotting results.")
print("Value sums:")
# TODO: come up with new, better plotting that puts them all in one plot instead of doing 1v1-comparison
for a1 in range(len(algorithm_list)):
    print(algorithm_list[a1].name, ": ", [np.round(z, 3) for z in algorithm_list[a1].values_sums])
    for a2 in range(a1 + 1, len(algorithm_list)):
        plotting.get_barchart_win_lose_draw(algo1_name=algorithm_list[a1].name, algo2_name=algorithm_list[a2].name, results_dict=results_dict, seed=seed, note="qbsolv")
print("Done plotting results.")
