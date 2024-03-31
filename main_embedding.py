import numpy as np
import random
import pickle
import time
import datetime
import os

import subprocess

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
    utils.append_to_pickle(data, f"data/tests/data_{graph_sizes}_{num_graphs_per_size}_{seed}.pkl")
    return data, True


def run_embedding_algorithm(serialized_algorithm, serialized_edges, num_agents, serialized_run_id, serialized_directory_path, graph_num):
    try:
        algorithm = pickle.loads(serialized_algorithm)
        edges = pickle.loads(serialized_edges)
        print(f"                    Measuring embedding of {algorithm.name} for graph_size {num_agents} for graph {graph_num} ...")
        algorithm.measure_embedding_run(num_agents, edges)
    except Exception as e:
        # for e.g. "No embedding found" exception
        raise Exception(f"{str(e)}") from None


def main(algorithm_list, data, graph_sizes, num_graphs_per_size, experiment, directory):
    run_id = str(datetime.datetime.now().date()) + '_' + str(datetime.datetime.now().time()).replace(':', '-')

    for algorithm in algorithm_list:
        too_large = False
        if isinstance(algorithm, IterativeQuantumAlgorithmWithK):
            agents = graph_sizes
            for num_agents in agents:
                if algorithm.k <= num_agents:
                    print(f"\n\n\nTest for graphsize {num_agents}")
                    this_range = range(num_graphs_per_size)
                    for graph_num in this_range:
                        print(f"\n\n     Graph {graph_num}")
                        graph = data[num_agents][graph_num]
                        if synthetic:
                            edges = graph
                        else:
                            edges = utils.transform(graph)

                        print(f"\n          Running {algorithm.name}...")
                        try:
                            # Run algorithm in a subprocess to avoid a potential SIGKILL of the entire script
                            # due to too much memory usage of this algorithm
                            # Serialize the objects and strings to pass to sub-process
                            serialized_algorithm = pickle.dumps(algorithm)
                            serialized_graph = pickle.dumps(edges)
                            serialized_run_id = pickle.dumps(run_id)
                            serialized_directory = pickle.dumps(directory)
                            # subprocess command
                            command = [
                                "python",
                                "-c",
                                f"from main_embedding import run_embedding_algorithm; "
                                f"print(run_embedding_algorithm({serialized_algorithm}, {serialized_graph}, {num_agents}, {serialized_run_id}, {serialized_directory}, {graph_num}))"
                            ]
                            # Run the subprocess with the memory limit
                            subprocess.run(command, check=True)
                        except subprocess.CalledProcessError as e:
                            print("                    Error: No embedding found, message: ", e)
                            too_large = True
                            break
                if too_large:
                    break
                else:
                    pass
        else:
            agents = graph_sizes
            for num_agents in agents:
                print(f"\n\n\nTest for graphsize {num_agents}")
                this_range = range(num_graphs_per_size)
                for graph_num in this_range:
                    print(f"\n\n     Graph {graph_num}")
                    graph = data[num_agents][graph_num]
                    if synthetic:
                        edges = graph
                    else:
                        edges = utils.transform(graph)

                    print(f"\n          Running {algorithm.name}...")
                    try:
                        # Run algorithm in a subprocess to avoid a potential SIGKILL of the entire script
                        # due to too much memory usage of this algorithm
                        # Serialize the objects and strings to pass to sub-process
                        serialized_algorithm = pickle.dumps(algorithm)
                        serialized_graph = pickle.dumps(edges)
                        serialized_run_id = pickle.dumps(run_id)
                        serialized_directory = pickle.dumps(directory)
                        # subprocess command
                        command = [
                            "python",
                            "-c",
                            f"from main_embedding import run_embedding_algorithm; "
                            f"print(run_embedding_algorithm({serialized_algorithm}, {serialized_graph}, {num_agents}, {serialized_run_id}, {serialized_directory}, {graph_num}))"
                        ]
                        # Run the subprocess with the memory limit
                        subprocess.run(command, check=True)
                    except subprocess.CalledProcessError as e:
                        print("                    Error: No embedding found, message: ", e)
                        too_large = True
                        break
                if too_large:
                    break
                else:
                    pass


if __name__ == "__main__":
    # TODO: Determine sensible number of seeds for statistical significance
    num_seeds = 1
    for _ in range(num_seeds):
        # Setting the seed
        #seed = random.randint(0, 2 ** 32 - 1)
        seed = 0  # Seed not relevant for D-Wave
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

        # D-Wave -> uncomment this and comment simulate for running with D-Wave
        solvers = ["dwave"]
        parallel = [True]
        k_list = [4, 3, 5, 2]

        for solver in solvers:
            '''
            algorithm_list = [#Jonas(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver),
                              Danielle(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver),
                              n_split_GCSQ(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver),
                              r_qubo(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver)
                              ]

            directory = f"results/{data_name}/quantum/{solver}"
            directory_exists = create_nested_directory(directory)
            if directory_exists:
                main(algorithm_list=algorithm_list, data=data, graph_sizes=graph_sizes, num_graphs_per_size=num_graphs_per_size,
                     experiment=f"non-iterative algorithms with {solver}", directory=directory)
            '''
            for mode in parallel:
                algorithm_list = [GCSQ(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode)]
                directory = f"results/{data_name}/quantum/{solver}/{'parallel' if mode else 'sequential'}"
                directory_exists = create_nested_directory(directory)
                if directory_exists:
                    main(algorithm_list=algorithm_list, data=data, graph_sizes=graph_sizes,
                         num_graphs_per_size=num_graphs_per_size, experiment=f"GCS-Q with {solver} in {mode} mode",
                         directory=directory)
                for k in k_list:
                    algorithm_list = [ours_iterative_exactly(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      ours_iterative_at_most(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      k_split_GCSQ_exactly(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      k_split_GCSQ_at_most(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k),
                                      r_qubo_iterative(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, parallel=mode, k=k)
                                      ]
                    directory = f"results/{data_name}/quantum/{solver}/{'parallel' if mode else 'sequential'}/k={k}"
                    directory_exists = create_nested_directory(directory)
                    if directory_exists:
                        main(algorithm_list=algorithm_list, data=data, graph_sizes=graph_sizes,
                             num_graphs_per_size=num_graphs_per_size, experiment=f"k-split algorithms with {solver} in {mode} mode for k={k}",
                             directory=directory)