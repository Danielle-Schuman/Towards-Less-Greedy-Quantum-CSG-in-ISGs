import utils

from algorithm import Algorithm


class GCSQ(Algorithm):
    def __init__(self, seed, num_graph_sizes, solver="qbsolv", timeout=10):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, timeout=timeout)
        self.name = f"GCS-Q_{self.solver}"

    def _split(self, c, edges):
        Q = {}
        for i in range(len(c)):
            for j in range(len(c)):
                if i < j:
                    utils.add(Q, i, i, edges[(c[i],c[j])])
                    utils.add(Q, j, j, edges[(c[i],c[j])])
                    utils.add(Q, i, j, -2*edges[(c[i],c[j])])
        solution = self.solve_qubo(Q, len(c))
        c1 = [c[k] for k in range(len(c)) if solution[k] == 1]
        c2 = [c[k] for k in range(len(c)) if solution[k] == 0]
        return [c1, c2]


