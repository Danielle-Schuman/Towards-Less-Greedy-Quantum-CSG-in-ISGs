import numpy as np
import random

import utils
import plotting
from algorithm import Algorithm
from GCSQ import GCSQ
from jonas import Jonas
from danielle import Danielle

# Setting the seed
seed = random.randint(0, 2**32 - 1)
random.seed(seed)
np.random.seed(seed)
print(f"Seed: {seed}")

tests = [10, 15, 20]
algorithm_list = [GCSQ(seed=seed), Jonas(seed=seed, num_coalitions=6), Danielle(seed=seed)]
results_dict = {}

for n in tests:
    print(f"Test for graphsize {n}")
    results_dict[n] = {}
    for graph_num in range(25):
        print(f"     Graph {graph_num}")
        edges = utils.generate_problem(n, mean=0.5)
        for algorithm in algorithm_list:
            coalitions = algorithm.solve(n, edges)
            v = np.sum([utils.value(c, edges) for c in coalitions])
            results_dict[n].setdefault(algorithm.name, []).append(v)
            print(f"          Coalition structure value for {algorithm.name}: {v}")
    print(f"Results for Test with graph size {n}: {results_dict[n]}\n")
print("Done running tests, plotting results.")
# TODO: come up with new, better plotting that puts them all in one plot instead of doing 1v1-comparison
for a1 in range(len(algorithm_list)):
    for a2 in range(a1 + 1, len(algorithm_list)):
        plotting.get_barchart_win_lose_draw(algo1_name=algorithm_list[a1].name, algo2_name=algorithm_list[a2].name, results_dict=results_dict, seed=seed, note="qbsolv")
print("Done plotting results.")
