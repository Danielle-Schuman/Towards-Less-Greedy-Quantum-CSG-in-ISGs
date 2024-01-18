import numpy as np
import random
import pickle
import time
import datetime
import os

import utils
from jonas import Jonas
from danielle import Danielle
from n_split_GCSQ import n_split_GCSQ
from r_qubo import r_qubo
from GCSQ import GCSQ
from algorithm import IterativeQuantumAlgorithmWithK
from ours_iterative_exactly import ours_iterative_exactly
from ours_iterative_at_most import ours_iterative_at_most
from k_split_GCSQ_exactly import k_split_GCSQ_exactly
from k_split_GCSQ_at_most import k_split_GCSQ_at_most
from r_qubo_iterative import r_qubo_iterative


def create_nested_directory(path):
    try:
        # Check if the path exists
        if not os.path.exists(path):
            # If not, create the missing directories
            os.makedirs(path)
            print(f"Created directory: {path}")
        return True
    except Exception as e:
        print(f"Error creating directory: {e}")
        return False

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


def main(algorithm_list, data, graph_sizes, num_graphs_per_size, experiment, directory, seed):
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
                if isinstance(algorithm, IterativeQuantumAlgorithmWithK):
                    if algorithm.k <= num_agents:
                        print(f"\n          Running {algorithm.name}...")
                        start_time = time.time()
                        coalitions = algorithm.solve(num_agents, edges)
                        end_time = time.time()
                        value = np.sum([utils.value(c, edges) for c in coalitions])
                        total_time = end_time - start_time
                        algorithm.data.append((coalitions, value, total_time))
                        pickle.dump(algorithm.data, open(f"{directory}/data_{algorithm.name}__{seed}__{run_id}.pkl", 'wb'))
                        results_dict[num_agents].setdefault(algorithm.name, []).append(value)
                        print(f"          Coalition structure value for {algorithm.name}: {value}   -   Time: {total_time}")
                        algorithm.values_sums[i] += value
                    else:
                        pass
                else:
                    print(f"\n          Running {algorithm.name}...")
                    start_time = time.time()
                    coalitions = algorithm.solve(num_agents, edges)
                    end_time = time.time()
                    value = np.sum([utils.value(c, edges) for c in coalitions])
                    total_time = end_time - start_time
                    algorithm.data.append((coalitions, value, total_time))
                    pickle.dump(algorithm.data, open(f"{directory}/data_{algorithm.name}_{seed}_{run_id}.pkl", 'wb'))
                    results_dict[num_agents].setdefault(algorithm.name, []).append(value)
                    print(f"          Coalition structure value for {algorithm.name}: {value}   -   Time: {total_time}")
                    algorithm.values_sums[i] += value
        print(f"Results for Test with graph size {num_agents}: {results_dict[num_agents]}\n")
    print(f"Done running tests for {experiment}.")

    '''
    # TODO: Move to plotting.py
    print("Value sums:")
    # TODO: come up with new, better plotting that puts them all in one plot instead of doing 1v1-comparison
    for a1 in range(len(algorithm_list)):
        print(algorithm_list[a1].name, ": ", [np.round(z, 3) for z in algorithm_list[a1].values_sums])
        for a2 in range(a1 + 1, len(algorithm_list)):
            plotting.get_barchart_win_lose_draw(algo1_name=algorithm_list[a1].name, algo2_name=algorithm_list[a2].name, results_dict=results_dict, seed=seed, note="qbsolv")
    print("Done plotting results.")
    '''


if __name__ == "__main__":
    # TODO: Determine sensible number of seeds for statistical significance
    num_seeds = 1
    for _ in range(num_seeds):
        # Setting the seed
        seed = random.randint(0, 2 ** 32 - 1)
        random.seed(seed)
        np.random.seed(seed)
        print(f"Seed: {seed}")

        # loading E.ON data
        data, synthetic = pickle.load(open('data/data_new_20samples_4_28.pkl', 'rb')), False
        # alternative: load synthetic data
        # data, synthetic = create_data_synthetic_test([10, 15, 20], 25, seed)
        if synthetic:
            data_name = "synthetic"
        else:
            data_name = "eon_data"

        graph_sizes = list(data.keys())  # is [4,6,8,10,12,14,16,18,20,22,24,26,28] for E.ON data
        num_graphs_per_size = len(data[graph_sizes[0]])  # is 20 for E.ON data
        num_graph_sizes = len(graph_sizes)

        # Simulate

        solvers = ["qbsolv", "sa", "qaoa"]
        parallel = [True, False]
        k_list = [i for i in range(2, graph_sizes[-1]+1)]

        # D-Wave -> uncomment this and comment simulate for running with D-Wave
        '''
        solvers = ["dwave"]
        parallel = [True]
        k_list = [4]  # TODO: Come up with sensible values based on simulation
        '''

        for solver in solvers:
            if solver == "qaoa":
                utils.qaoa_init(seed=seed)
            if solver == "dwave":
                utils.dwave_init()
            algorithm_list = [Jonas(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver)] #,
                              #Danielle(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver),
                              #n_split_GCSQ(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver),
                              #r_qubo(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver)]
            directory = f"results/{data_name}/quantum/{solver}"
            directory_exists = create_nested_directory(directory)
            if directory_exists:
                main(algorithm_list=algorithm_list, data=data, graph_sizes=graph_sizes, num_graphs_per_size=num_graphs_per_size,
                     experiment=f"non-iterative algorithms with {solver}", directory=directory, seed=seed)
            for mode in parallel:
                algorithm_list = [GCSQ(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode)]
                directory = f"results/{data_name}/quantum/{solver}/{'parallel' if mode else 'sequential'}"
                directory_exists = create_nested_directory(directory)
                if directory_exists:
                    main(algorithm_list=algorithm_list, data=data, graph_sizes=graph_sizes,
                         num_graphs_per_size=num_graphs_per_size, experiment=f"GCS-Q with {solver} in {mode} mode",
                         directory=directory, seed=seed)
                for k in k_list:
                    algorithm_list = [ours_iterative_exactly(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      ours_iterative_at_most(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      k_split_GCSQ_exactly(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      k_split_GCSQ_at_most(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      r_qubo_iterative(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k)]
                    directory = f"results/{data_name}/quantum/{solver}/{'parallel' if mode else 'sequential'}/k={k}"
                    directory_exists = create_nested_directory(directory)
                    if directory_exists:
                        main(algorithm_list=algorithm_list, data=data, graph_sizes=graph_sizes,
                             num_graphs_per_size=num_graphs_per_size, experiment=f"k-split algorithms with {solver} in {mode} mode for k={k}",
                             directory=directory, seed=seed)

        # Classical algorithms
        '''
        algorithm_list = []  # TODO
        directory = f"results/{data_name}/classical"
        directory_exists = create_nested_directory(directory)
        if directory_exists:
            main(algorithm_list=algorithm_list, data=data, graph_sizes=graph_sizes,
                 num_graphs_per_size=num_graphs_per_size,
                 experiment=f"classical algorithms",
                 directory=directory, seed=seed)
        '''

