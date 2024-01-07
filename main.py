import numpy as np
import random
import pickle
import time

import utils
import plotting
from algorithm import Algorithm
from GCSQ import GCSQ
from jonas import Jonas
from danielle import Danielle
from min_k_cut_non_iterative import min_k_cut_non_iterative
from ours_iterative_exactly import ours_iterative_exactly
from ours_iterative_at_most import ours_iterative_at_most
from min_k_cut_exactly import min_k_cut_exactly
from min_k_cut_at_most import min_k_cut_at_most

# Setting the seed
seed = random.randint(0, 2**32 - 1)
random.seed(seed)
np.random.seed(seed)
print(f"Seed: {seed}")

# loading data
#data = pickle.load(open('data_new_20samples_4_28.pkl', 'rb'))

tests = [10, 15, 20]
# TODO: Find sensible k's to try, possibly everything between 2 and n
k = 4
algorithm_list = [GCSQ(seed=seed), Jonas(seed=seed, num_coalitions=6), Danielle(seed=seed), min_k_cut_non_iterative(seed=seed),
                  ours_iterative_exactly(seed=seed, k=k), ours_iterative_at_most(seed=seed, k=k), min_k_cut_exactly(seed=seed, k=k), min_k_cut_at_most(seed=seed, k=k)]
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
