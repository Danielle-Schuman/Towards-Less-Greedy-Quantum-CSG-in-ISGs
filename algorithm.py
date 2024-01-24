import copy
from abc import ABC, abstractmethod
import utils
import pickle


class Algorithm(ABC):
    def __init__(self, seed, num_graph_sizes, category, num_coalitions=None, timeout=2592000):
        self.category = category
        self.num_coalitions = num_coalitions
        self.seed = seed
        self.timeout = timeout
        # list of tupel (coalitions, value, total_time) for last graph evaluated
        # later for plotting: list of all such tuples
        self.data = None
        # list of sums of values for graphs of same size
        self.values_sums = [0] * num_graph_sizes

        self.name = "missing_algorithm_name"  # gets overwritten with the actual name in each of the classes

    @abstractmethod
    def solve(self, num_agents, edges):
        pass


class ClassicalAlgorithm(Algorithm, ABC):
    def __init__(self, seed, num_graph_sizes, num_coalitions=None, timeout=2592000):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, category="classical", num_coalitions=num_coalitions,
                         timeout=timeout)


class QuantumAlgorithm(Algorithm, ABC):
    def __init__(self, seed, num_graph_sizes, solver="qbsolv", num_coalitions=None, timeout=2592000):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, category="quantum", num_coalitions=num_coalitions,
                         timeout=timeout)
        self.solver = solver  # can be qbsolv, sa, dwave or qaoa

    def solve_qubo(self, qubo, num_qubits):
        if self.solver == "qbsolv":
            solution = utils.solve_with_qbsolv(qubo, num_qubits, self.seed, self.timeout)
        elif self.solver == "sa":
            solution = utils.solve_with_sa(qubo, num_qubits, self.seed)
        elif self.solver == "qaoa":
            solution = utils.solve_with_qaoa(qubo, num_qubits, self.seed)
        elif self.solver == "dwave":
            solution = utils.solve_with_dwave(qubo, num_qubits, self.solver)
        else:
            print("Chosen solver not available. Solving with QB-Solv.")
            solution = utils.solve_with_qbsolv(qubo, num_qubits, self.seed, self.timeout)
        return solution


class IterativeQuantumAlgorithm(QuantumAlgorithm, ABC):
    def __init__(self, seed, num_graph_sizes, solver="qbsolv", num_coalitions=None, timeout=2592000, parallel=True):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, num_coalitions=num_coalitions, timeout=timeout)
        self.parallel = parallel

    @staticmethod
    def convert_qubo(qubo, old_num_qubits):
        updated_qubo = {}
        for key, value in qubo.items():
            updated_key = tuple(num + old_num_qubits for num in key)
            updated_qubo[updated_key] = value
        return updated_qubo

    @staticmethod
    def split_solution(solution, qubo_starts):
        result_lists = []
        # Add sublists between consecutive qubo_starts
        for i in range(len(qubo_starts) - 1):
            start_index = qubo_starts[i]
            end_index = qubo_starts[i + 1]
            result_lists.append(solution[start_index:end_index])
        # Add the last sublist from the last qubo_start to the end
        result_lists.append(solution[qubo_starts[-1]:])
        return result_lists

    def solve_qubos(self, qubo_list, coalitions):
        converted_qubo = qubo_list[0][0]
        total_num_qubits = qubo_list[0][1]
        qubit_num_start_of_qubo = [0]
        for (qubo, num_qubits) in qubo_list[1:]:
            qubit_num_start_of_qubo.append(total_num_qubits)
            qubo_with_new_numbers = IterativeQuantumAlgorithm.convert_qubo(qubo, total_num_qubits)
            total_num_qubits = total_num_qubits + num_qubits
            converted_qubo.update(qubo_with_new_numbers)
        all_solutions = self.solve_qubo(converted_qubo, total_num_qubits)
        solutions = IterativeQuantumAlgorithm.split_solution(all_solutions, qubit_num_start_of_qubo)
        new_coalitions = []
        for i, solution in enumerate(solutions):
            new_coalitions.append(self._get_coalitions_from_qubo_solution(coalitions[i], solution))
        return new_coalitions

    def _split(self, coalition, edges):
        Q, num_qubits = self._get_qubo(coalition, edges)
        solution = self.solve_qubo(Q, num_qubits)
        return self._get_coalitions_from_qubo_solution(coalition, solution)

    # implementation is a tiny bit different from the one in the GCS-Q paper, but does the same thing
    def solve(self, num_agents, edges):
        # initialize with grand coalition
        coalitions = [list(range(num_agents))]
        # for a maximum of num_agents steps (as we can at most have num_agents coalitions, one for each agent)
        for n in range(num_agents):
            new_coalitions = copy.deepcopy(coalitions)
            # exclude coalitions that cannot be further split because they already contain only one agent
            real_coalitions = [coalition for coalition in coalitions if len(coalition) > 1]
            if self.parallel:
                qubo_list = []
            for coalition in real_coalitions:
                if self.parallel:
                    qubo, num_qubits = self._get_qubo(coalition, edges)
                    qubo_list.append((qubo, num_qubits))
                else:
                    # Solve optimal split problem
                    # Let annealer find optimal split of coalition
                    split_result = self._split(coalition, edges)
                    # If splitting increases value of this part of the coalition structure,
                    # replace coalition with split coalitions
                    if sum([utils.value(c_new, edges) for c_new in split_result]) > utils.value(coalition, edges):
                        new_coalitions.remove(coalition)
                        new_coalitions.extend(split_result)
            if self.parallel:
                split_results = self.solve_qubos(qubo_list, real_coalitions)
                for i, split_result in enumerate(split_results):
                    # If splitting increases value of this part of the coalition structure,
                    # replace coalition with split coalitions
                    if sum([utils.value(c_new, edges) for c_new in split_result]) > utils.value(real_coalitions[i], edges):
                        new_coalitions.remove(real_coalitions[i])
                        new_coalitions.extend(split_result)
            # none of the splittings of the coalitions of the previous coalition structure improved it, stop the algorithm
            # and return the previous coalition structure
            if len(coalitions) == len(new_coalitions):
                break
            coalitions = new_coalitions
            if self.solver == "dwave":
                pickle.dump(coalitions, open(f"results/coalitions_{self.name}_after_n={n}.pkl", 'wb'))
        return coalitions

    @abstractmethod
    def _get_qubo(self, coalition, edges):
        pass

    @abstractmethod
    def _get_coalitions_from_qubo_solution(self, coalition, solution):
        pass


class IterativeQuantumAlgorithmWithK(IterativeQuantumAlgorithm, ABC):
    def __init__(self, seed, num_graph_sizes, k, solver="qbsolv", num_coalitions=None, timeout=2592000, parallel=True):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, num_coalitions=num_coalitions, timeout=timeout, parallel=parallel)
        self.k = k
