import copy
import numpy as np
import pickle
from dwave_qbsolv import QBSolv
import dimod as di
from neal import SimulatedAnnealingSampler
from uqo.client.config import Config
from uqo import Problem
from qiskit_optimization import QuadraticProgram
from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_algorithms.utils import algorithm_globals
from qiskit.primitives import Sampler


# this function solves a given QUBO-Matrix qubo with QB-solv
# see: https://docs.ocean.dwavesys.com/projects/qbsolv/en/latest/source/generated/dwave_qbsolv.QBSolv.sample.html
def solve_with_qbsolv(qubo, num_qubits, seed, timeout=2592000):
    # num_repeats: Determines the number of times to repeat the main loop in qbsolv after determining a better sample. Default 50.
    # timeout: Number of seconds before routine halts. Default is 2592000 = 1 month. Jonas has timeout=10, but that can marginally decrease performance.
    # actually, this has sample_count = 1 -> for sample_count = 100, we'd need to run it 100 times with different seeds
    # but it seems to perform just fine without
    response = QBSolv().sample_qubo(qubo, num_repeats=1000, timeout=timeout, seed=seed)
    solution = [response.samples()[0][i] for i in range(num_qubits)]
    return solution


# see: https://docs.ocean.dwavesys.com/projects/neal/en/latest/reference/generated/neal.sampler.SimulatedAnnealingSampler.sample.html
# (timeout with interrupt_function seems to degenerate performance, even when not hitting the timeout, for some reason -> not used)
def solve_with_sa(qubo, num_qubits, seed):
    shots = 100  # number of times the SA algorithms is executed, default is 1
    anneal_steps = 1000  # number of updates of "qubits values" in SA algorithm, default is 1000
    qubo_as_bqm = di.BQM(qubo, "BINARY")
    response = SimulatedAnnealingSampler().sample(qubo_as_bqm, num_reads=shots, num_sweeps=anneal_steps, seed=seed)
    solution = [response.samples()[0][i] for i in range(num_qubits)]
    return solution


def solve_with_dwave(qubo, num_qubits, solver):
    # get advantage solver for connection
    config = Config(configpath="secrets/config.json")
    connection = config.create_connection()
    available_solvers = connection.get_available_dwave_solvers()
    if 'Advantage_system6.3' in available_solvers:
        advantage_solver = 'Advantage_system6.3'
    elif 'Advantage_system4.1' in available_solvers:
        advantage_solver = 'Advantage_system4.1'
    else:
        raise Exception("No know Advantage solver available.")
    print(f"Running on {advantage_solver} ...")

    shots = 100  # number of times the QA algorithms is executed, default is 1 -> TODO: Find good value
    if solver == "test_uqo":
        # qbsolv for testing uqo connection without decreasing our quota
        response = Problem.Qubo(config, qubo).with_platform("qbsolv").solve(shots)
    elif solver == "dwave":
        global ADVANTAGE_SOLVER
        # needs dwave quota
        problem = Problem.Qubo(config, qubo).with_platform("dwave").with_solver(advantage_solver)
        # calculate embedding
        problem.find_pegasus_embedding()
        try:
            response = problem.solve(shots)
        except:
            response = problem.solve(shots)
    solution = [response.solutions[0][i] for i in range(num_qubits)]
    return solution


def convert_dict_keys_to_strings(input_dict):
    new_dict = {}
    for key, value in input_dict.items():
        new_key = tuple(str(i) for i in key)
        new_dict[new_key] = value
    return new_dict

def convert_dict_to_list(input_dict):
    max_key = max(map(int, input_dict.keys())) if input_dict else 0
    result_list = [None] * (max_key + 1)
    for key, value in input_dict.items():
        result_list[int(key)] = value
    return result_list


def solve_with_qaoa(qubo, num_qubits, seed):
    # init qaoa algorithm
    # 100 shots would be "fair", but doesn't find the optimum even for small examples
    shots = 1000
    algorithm_globals.random_seed = seed
    # classical Optimizer (uses Machine Learning)
    # maxiter = num circuit evaluations -> a bit like num_repeats in QBsolv or num_sweeps in SA, default=1000
    # (100 works too, but not much faster than 1000; 10 degrades performance)
    optimizer = COBYLA(maxiter=1000, disp=False)  # disp = Whether it prints out information
    # has "None" as default for shots -> in that case, it calculates probabilities
    sampler = Sampler()
    sampler.set_options(shots=shots, seed=seed)
    qaoa = QAOA(optimizer=optimizer, sampler=sampler)
    qaoa_algorithm = MinimumEigenOptimizer(qaoa)

    # Convert qubo to right format
    qubo_for_qaoa = QuadraticProgram('qubo')
    for n in range(num_qubits):
        qubo_for_qaoa.binary_var(name=str(n))
    qubo = convert_dict_keys_to_strings(qubo)
    qubo_for_qaoa.minimize(quadratic=qubo)
    # Run quantum algorithm QAOA
    result = qaoa_algorithm.solve(qubo_for_qaoa)
    return convert_dict_to_list(result.variables_dict)


def append_to_pickle(data, file_path):
    try:
        with open(file_path, 'ab') as file:
            pickle.dump(data, file)
    except FileNotFoundError:
        with open(file_path, 'wb') as file:
            pickle.dump(data, file)
        print(f"Created {file_path} and wrote first data.")
    except Exception as e:
        print(f"An error occurred trying to pickle the data {data} in the file {file_path}:", e)


# this function calculates the value of a solution for a given QUBO-Matrix qubo
def getValue(Q, solution):
    ones = [i for i in range(len(solution)) if solution[i] == 1]
    value = 0
    for x in ones:
        for y in ones:
            if (x,y) in Q.keys():
                value += Q[(x,y)]
    return value


def transform(naeimeh_graph):
    edges = {}
    for key in naeimeh_graph.keys():
        i, j = key.split(",")
        i, j = int(i)-1, int(j)-1
        edges[(i,j)] = naeimeh_graph[key]
    return edges


# this function prints the first num_agents row/columns of a QUBO-Matrix qubo
def printQUBO(Q, n):
    for row in range(n):
        for column in range(n):
            if row > column:
                print("        ", end = '')
                continue
            printing = ""
            if (row,column) in Q.keys() and Q[(row,column)] != 0:
                printing = str(float(Q[(row,column)]))
            printing += "_______"
            printing = printing[:7]
            printing += " "
            print(printing, end = '')
        print("")



def generate_problem(n, mean=0.5):
    edges = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                edges[(i,j)] = np.random.sample() - mean
    return edges



def add(Q, i, j, v):
    if (i,j) not in Q.keys():
        Q[(i,j)] = v
    else:
        Q[(i,j)] += v



def value(c, edges):
    v = 0
    for i in c:
        for j in c:
            if i < j:
                v += edges[(i,j)]
    return v


def baseline(n, edges):

    best_coalitions = None
    best_value = 0

    for _ in range(10):
        agent_list = np.array(range(n))
        np.random.shuffle(agent_list)
        coalitions = sub_baseline(edges, agent_list)
        v = np.sum([value(c,edges) for c in coalitions])
        if v > best_value:
            best_coalitions = coalitions
            best_value = v

    return best_coalitions


def sub_baseline(edges, agent_list):
    coalitions = [[agent_list[0]]]
    for i in agent_list[1:]:
        improvements = []
        for c in coalitions:
            old_v = value(c, edges)
            new_c = copy.deepcopy(c)
            new_c.append(i)
            new_v = value(new_c, edges)
            improvements.append(new_v - old_v)
        if np.max(improvements) > 0:
            coalitions[np.argmax(improvements)].append(i)
        else:
            coalitions.append([i])
    return coalitions