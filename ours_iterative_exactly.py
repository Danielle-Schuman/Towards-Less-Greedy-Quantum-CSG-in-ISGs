import utils
import numpy as np

from algorithm import IterativeAlgorithmWithK


class ours_iterative_exactly(IterativeAlgorithmWithK):
    def __init__(self, seed, num_graph_sizes, k, solver="qbsolv", timeout=10, parallel=True):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, timeout=timeout, k=k, parallel=True)
        self.name = f"{self.k}_split_ours_iterative_exactly_{self.solver}_{'parallel' if self.parallel else 'sequential'}"

    def _get_qubo(self, coalition, edges):
        Q = {}
        # TODO: Find good value for penalty using penalty engineering
        # sum of all the edges absolute values to use as a penalty value later
        penalty = np.sum(np.abs(list(edges.values())))
        # iterate over agents vertically (rows)
        for i in range(len(coalition)):
            for c in range(self.k):
                # get number of logical qubit vertically (rows)
                q_ic = i * self.k + c
                # iterate over agents horizontally (columns)
                for j in range(i + 1, len(coalition)):
                    # get number of logical qubit horizontally (columns)
                    q_jc = j * self.k + c
                    # put the negative of the edge weight in the graph as an incentive
                    # to put the agents in the same coalition if the edge weight is > 0
                    utils.add(Q, q_ic, q_jc, -edges[(coalition[i], coalition[j])])
                # add reward for putting agent in any coalition
                utils.add(Q, q_ic, q_ic, -penalty)
                for c2 in range(c + 1, self.k):
                    # get number of logical qubit horizontally (columns)
                    q_ic2 = i * self.k + c2
                    # add penalty for putting agent in two different coalitions at the same time (we don't do overlap here yet)
                    utils.add(Q, q_ic, q_ic2, 2 * penalty)
        return Q, self.k * len(coalition)

    def _get_coalitions_from_qubo_solution(self, coalition, solution):
        # make a list of coalitions, with each coalition being a list with the numbers of the agents in these coalitions
        coalitions = []
        for c in range(self.k):
            new_coalition = [coalition[i] for i in range(len(coalition)) if solution[i * self.k + c] == 1]
            coalitions.append(new_coalition)
        return coalitions
    