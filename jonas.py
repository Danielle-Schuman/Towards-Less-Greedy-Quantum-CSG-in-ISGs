import utils
import numpy as np

from algorithm import Algorithm


class Jonas(Algorithm):
    def __init__(self, seed,  num_graph_sizes, solver="qbsolv", num_coalitions=None, timeout=10):
        super().__init__(seed=seed,  num_graph_sizes=num_graph_sizes, solver=solver, num_coalitions=num_coalitions, timeout=timeout)
        self.name = f"ours_n_half_{self.solver}"

    def solve(self, num_agents, edges):
        if not self.num_coalitions:  # self.num_coalitions is still None
            self.num_coalitions = num_agents // 2
        # create the QUBO
        Q = {}
        # sum of all the edges absolute values to use as a penalty value later
        penalty = np.sum(np.abs(list(edges.values())))
        # iterate over agents vertically (rows)
        for i in range(num_agents):
            # iterate over coalition groups vertically (rows)
            for c_i in range(self.num_coalitions):
                # get number of logical qubit vertically (rows)
                p_i = c_i * num_agents + i
                # iterate over agents horizontally (columns)
                for j in range(num_agents):
                    # iterate over coalition groups horizontally (columns)
                    for c_j in range(self.num_coalitions):
                        # get number of logical qubit horizontally (columns)
                        p_j = c_j * num_agents + j
                        # if we are in the upper triangular matrix looking at
                        # two different agents and the same coalition
                        if p_i < p_j and i != j and c_i == c_j:
                            # put the negative of the edge weight in the graph as an incentive
                            # to put the agents in the same coalition if the edge weight is > 0
                            utils.add(Q, p_i, p_j, -edges[(i,j)])
                        # if we are in the upper triangular matrix looking at
                        # the same agent and two different coalitions
                        elif p_i < p_j and i == j and c_i != c_j:
                            # add sum of all the edges absolute values as a penalty value to ensure that one agent cannot be in two coalitions at once
                            utils.add(Q, p_i, p_j, penalty)
        # solve the QUBO
        solution = self.solve_qubo(Q, self.num_coalitions * num_agents)
        # make a list of coalitions, with each coalition being a list with the numbers of the agents in these coalitions
        coalitions = []
        for c in range(self.num_coalitions):
            coalition = [k for k in range(num_agents) if solution[c * num_agents:(c + 1) * num_agents][k] == 1]
            coalitions.append(coalition)
        return coalitions

