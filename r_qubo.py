import utils
import numpy as np

from algorithm import QuantumAlgorithm


class r_qubo(QuantumAlgorithm):
    def __init__(self, seed, num_graph_sizes, solver="qbsolv", num_coalitions=None, timeout=2592000):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, num_coalitions=num_coalitions, solver=solver, timeout=timeout)
        self.name = f"R-QUBO_{self.solver}"

    def solve(self, num_agents, edges):
        if not self.num_coalitions:  # self.num_coalitions is still None
            self.num_coalitions = num_agents
        Q = {}
        # TODO: Find good value for penalty using penalty engineering
        # sum of all the edges absolute values to use as a penalty value later
        penalty = np.sum(np.abs(list(edges.values())))
        # iterate over agents vertically (rows)
        for i in range(num_agents):
            for c in range(self.num_coalitions - 1):
                # get number of logical qubit vertically (rows)
                q_ic = i * (self.num_coalitions - 1) + c
                # iterate over agents horizontally (columns)
                for j in range(i + 1, num_agents):
                    # get number of logical qubit horizontally (columns)
                    q_jc = j * (self.num_coalitions - 1) + c
                    utils.add(Q, q_ic, q_ic, edges[(i, j)])
                    utils.add(Q, q_jc, q_jc, edges[(i, j)])
                    utils.add(Q, q_ic, q_jc, -2 * edges[(i, j)])
                    for c2 in range(self.num_coalitions - 1):
                        if c2 != c:
                            # get number of logical qubit horizontally (columns)
                            q_jc2 = j * (self.num_coalitions - 1) + c2
                            utils.add(Q, q_ic, q_jc2, -edges[(i, j)])
            # put penalties for more than one coalition per agent
            for c1 in range(self.num_coalitions - 2):
                # get number of logical qubit vertically (rows)
                q_ic1 = i * (self.num_coalitions - 1) + c1
                for c2 in range(c1 + 1, (self.num_coalitions - 1)):
                    # get number of logical qubit horizontally (columns)
                    q_ic2 = i * (self.num_coalitions - 1) + c2
                    # add penalty for putting agent in any two coalitions
                    utils.add(Q, q_ic1, q_ic2, penalty)

        # solve the QUBO
        solution = self.solve_qubo(Q, (self.num_coalitions - 1) * num_agents)
        # make a list of coalitions, with each coalition being a list with the numbers of the agents in these coalitions
        coalitions = []
        for c in range(self.num_coalitions - 1):
            coalition = [i for i in range(num_agents) if solution[i * (self.num_coalitions - 1) + c] == 1]
            coalitions.append(coalition)
        coalition = [i for i in range(num_agents) if solution[(i * (self.num_coalitions - 1)):((i + 1) * (self.num_coalitions - 1))] == ([0] * (self.num_coalitions - 1))]
        coalitions.append(coalition)
        return coalitions

