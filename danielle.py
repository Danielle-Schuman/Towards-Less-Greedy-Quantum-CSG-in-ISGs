import utils
import numpy as np

from algorithm import QuantumAlgorithm


class Danielle(QuantumAlgorithm):
    def __init__(self, seed, num_graph_sizes, solver="qbsolv", num_coalitions=None, timeout=2592000):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, solver=solver, num_coalitions=num_coalitions, timeout=timeout)
        self.name = f"ours_n_{self.solver}"

    def solve(self, num_agents, edges):
        if not self.num_coalitions:  # self.num_coalitions is still None
            self.num_coalitions = num_agents
        # create the QUBO
        Q = {}
        # TODO: Find good value for penalty using penalty engineering
        # sum of all the edges absolute values to use as a penalty value later
        penalty = np.sum(np.abs(list(edges.values())))
        # iterate over agents vertically (rows)
        for i in range(num_agents):
            for c in range(self.num_coalitions):
                # get number of logical qubit vertically (rows)
                q_ic = i * self.num_coalitions + c
                # iterate over agents horizontally (columns)
                for j in range(i + 1, num_agents):
                    # get number of logical qubit horizontally (columns)
                    q_jc = j * self.num_coalitions + c
                    # put the negative of the edge weight in the graph as an incentive
                    # to put the agents in the same coalition if the edge weight is > 0
                    utils.add(Q, q_ic, q_jc, -edges[(i, j)])
                # add reward for putting agent in any coalition
                utils.add(Q, q_ic, q_ic, -penalty)
                for c2 in range(c + 1, self.num_coalitions):
                    # get number of logical qubit horizontally (columns)
                    q_ic2 = i * self.num_coalitions + c2
                    # add penalty for putting agent in two different coalitions at the same time (we don't do overlap here yet)
                    utils.add(Q, q_ic, q_ic2, 2 * penalty)

        # solve the QUBO
        solution = self.solve_qubo(Q, self.num_coalitions * num_agents)
        # make a list of coalitions, with each coalition being a list with the numbers of the agents in these coalitions
        coalitions = []
        for c in range(self.num_coalitions):
            coalition = [i for i in range(num_agents) if solution[i * self.num_coalitions + c] == 1]
            coalitions.append(coalition)
        return coalitions

