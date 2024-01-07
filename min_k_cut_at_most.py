import utils
import numpy as np

from algorithm import Algorithm


class min_k_cut_at_most(Algorithm):
    def __init__(self, seed, k, timeout=10):
        super().__init__(seed, timeout)
        self.name = "min_k_cut_at_most"
        self.k = k

    def _split(self, coalition, edges):
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
                    utils.add(Q, q_ic, q_ic, edges[(coalition[i], coalition[j])])
                    utils.add(Q, q_jc, q_jc, edges[(coalition[i], coalition[j])])
                    utils.add(Q, q_ic, q_jc, -2 * edges[(coalition[i], coalition[j])])
            # put penalties for more than one coalition per agent
            for c1 in range(self.k - 1):
                # get number of logical qubit vertically (rows)
                q_ic1 = i * self.k + c1
                for c2 in range(c1 + 1, self.k):
                    # get number of logical qubit horizontally (columns)
                    q_ic2 = i * self.k + c2
                    # add penalty for putting agent in any coalition
                    utils.add(Q, q_ic1, q_ic2, penalty)

        # solve the QUBO
        solution = utils.solve_with_qbsolv(Q, self.k * len(coalition), self.seed, self.timeout)
        # make a list of coalitions, with each coalition being a list with the numbers of the agents in these coalitions
        coalitions = []
        for c in range(self.k):
            new_coalition = [coalition[i] for i in range(len(coalition)) if solution[i * self.k + c] == 1]
            coalitions.append(new_coalition)
        return coalitions

