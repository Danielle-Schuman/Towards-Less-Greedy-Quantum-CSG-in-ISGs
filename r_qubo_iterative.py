import utils
import numpy as np

from algorithm import IterativeQuantumAlgorithmWithK


class r_qubo_iterative(IterativeQuantumAlgorithmWithK):
    def __init__(self, seed, num_graph_sizes, k, solver="qbsolv", timeout=600, parallel=True):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, timeout=timeout, k=k, parallel=parallel)
        self.name = f"{self.k}_split_R_QUBO-iterative_{self.solver}_{'parallel' if self.parallel else 'sequential'}"

    def _get_qubo(self, coalition, edges):
        Q = {}
        # TODO: Find good value for penalty using penalty engineering
        # sum of all the edges absolute values to use as a penalty value later
        penalty = np.sum(np.abs(list(edges.values())))
        # iterate over agents vertically (rows)
        for i in range(len(coalition)):
            for c in range(self.k - 1):
                # get number of logical qubit vertically (rows)
                q_ic = i * (self.k - 1) + c
                # iterate over agents horizontally (columns)
                for j in range(i + 1, len(coalition)):
                    # get number of logical qubit horizontally (columns)
                    q_jc = j * (self.k - 1) + c
                    utils.add(Q, q_ic, q_ic, edges[(coalition[i], coalition[j])])
                    utils.add(Q, q_jc, q_jc, edges[(coalition[i], coalition[j])])
                    utils.add(Q, q_ic, q_jc, -2 * edges[(coalition[i], coalition[j])])
                    for c2 in range(self.k - 1):
                        if c2 != c:
                            # get number of logical qubit horizontally (columns)
                            q_jc2 = j * (self.k - 1) + c2
                            utils.add(Q, q_ic, q_jc2, -edges[(coalition[i], coalition[j])])
            # put penalties for more than one coalition per agent
            for c1 in range(self.k - 2):
                # get number of logical qubit vertically (rows)
                q_ic1 = i * (self.k - 1) + c1
                for c2 in range(c1 + 1, (self.k - 1)):
                    # get number of logical qubit horizontally (columns)
                    q_ic2 = i * (self.k - 1) + c2
                    # add penalty for putting agent in any two coalitions
                    utils.add(Q, q_ic1, q_ic2, penalty)
        return Q, (self.k - 1) * len(coalition)

    def _get_coalitions_from_qubo_solution(self, coalition, solution):
        # make a list of coalitions, with each coalition being a list with the numbers of the agents in these coalitions
        coalitions = []
        for c in range(self.k - 1):
            new_coalition = [coalition[i] for i in range(len(coalition)) if solution[i * (self.k - 1) + c] == 1]
            coalitions.append(new_coalition)
        new_coalition = [coalition[i] for i in range(len(coalition)) if
                     solution[(i * (self.k - 1)):((i + 1) * (self.k - 1))] == (
                                 [0] * (self.k - 1))]
        coalitions.append(new_coalition)
        return coalitions

