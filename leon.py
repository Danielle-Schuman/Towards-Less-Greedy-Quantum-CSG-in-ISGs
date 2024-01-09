import utils

from algorithm import Algorithm


class Leon(Algorithm):
    def __init__(self, seed, num_graph_sizes, num_coalitions=None, timeout=10):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, num_coalitions=num_coalitions, timeout=timeout)
        self.name = "Leon19"

    def solve(self, num_agents, edges):
        #TODO
        pass
        