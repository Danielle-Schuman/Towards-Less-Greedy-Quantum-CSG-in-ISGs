import copy
import utils

class Algorithm:
    def __init__(self, seed, num_graph_sizes, num_coalitions=None, timeout=10):
        self.num_coalitions = num_coalitions
        self.seed = seed
        self.timeout = timeout
        # list of tupels (coalitions, value, total_time) for all graph evaluated so far
        self.data = []
        # list of sums of values for graphs of same size
        self.values_sums = [0] * num_graph_sizes


    # implementation is a tiny bit different from the one in the GCS-Q paper, but does the same thing
    def solve(self, num_agents, edges):
        # initialize with grand coalition
        coalitions = [list(range(num_agents))]
        # for a maximum of num_agents steps (as we can at most have num_agents coalitions, one for each agent)
        for _ in range(num_agents):
            new_coalitions = copy.deepcopy(coalitions)
            for c in coalitions:
                # unless the coalition cannot be further split because it already contains only one agent
                if len(c) > 1:
                    # Solve optimal split problem
                    # Let annealer find optimal split of coalition
                    split_result = self._split(c, edges)
                    # If splitting increases value of this part of the coalition structure,
                    # replace coalition c with split coalitions
                    if sum([utils.value(c_new, edges) for c_new in split_result]) > utils.value(c, edges):
                        new_coalitions.remove(c)
                        new_coalitions.extend(split_result)
            # none of the splittings of the coalitions of the previous coalition structure improved it, stop the algorithm
            # and return the previous coalition structure
            if len(coalitions) == len(new_coalitions):
                break
            coalitions = new_coalitions
        return coalitions