import utils

from algorithm import IterativeQuantumAlgorithm


class GCSQ(IterativeQuantumAlgorithm):
    def __init__(self, seed, num_graph_sizes, solver="qbsolv", timeout=10, parallel=True):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, timeout=timeout, parallel=parallel)
        self.name = f"GCS-Q_{self.solver}_{'parallel' if self.parallel else 'sequential'}"

    def _get_qubo(self, coalition, edges):
        Q = {}
        for i in range(len(coalition)):
            for j in range(len(coalition)):
                if i < j:
                    utils.add(Q, i, i, edges[(coalition[i], coalition[j])])
                    utils.add(Q, j, j, edges[(coalition[i], coalition[j])])
                    utils.add(Q, i, j, -2 * edges[(coalition[i], coalition[j])])
        return Q, len(coalition)

    def _get_coalitions_from_qubo_solution(self, coalition, solution):
        c1 = [coalition[k] for k in range(len(coalition)) if solution[k] == 1]
        c2 = [coalition[k] for k in range(len(coalition)) if solution[k] == 0]
        return [c1, c2]

