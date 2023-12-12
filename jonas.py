import utils

from algorithm import Algorithm


class Jonas(Algorithm):
    def __init__(self, seed, num_coalitions=3):
        super().__init__(seed, num_coalitions)
        self.name = "$Jonas_2$"

    def solve(self, n, edges):
        # create the QUBO
        Q = {}
        # iterate over agents vertically (rows)
        for i in range(n):
            # iterate over coalition groups vertically (rows)
            for c_i in range(self.num_coalitions):
                # get number of logical qubit vertically (rows)
                p_i = c_i * n + i
                # iterate over agents horizontally (columns)
                for j in range(n):
                    # iterate over coalition groups horizontally (columns)
                    for c_j in range(self.num_coalitions):
                        # get number of logical qubit horizontally (columns)
                        p_j = c_j * n + j
                        # if we are in the upper triangular matrix looking at
                        # two different agents and the same coalition
                        if p_i < p_j and i != j and c_i == c_j:
                            # put the negative of the edge weight in the graph as an incentive
                            # to put the agents in the same coalition if the edge weight is > 0
                            utils.add(Q, p_i, p_j, -edges[(i,j)])
                        # if we are in the upper triangular matrix looking at
                        # the same agent and two different coalitions
                        elif p_i < p_j and i == j and c_i != c_j:
                            # add 99999 as a penalty value to ensure that one agent cannot be in two coalitions at once
                            utils.add(Q, p_i, p_j, 99999)
        # solve the QUBO
        solution = utils.solve_with_qbsolv(Q, self.num_coalitions*n, self.seed)
        # make a list of coalitions, with each coalition being a list with the numbers of the agents in these coalitions
        coalitions = []
        for c in range(self.num_coalitions):
            coalition = [k for k in range(n) if solution[c*n:(c+1)*n][k] == 1]
            coalitions.append(coalition)
        return coalitions

