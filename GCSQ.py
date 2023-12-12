import utils
import copy

from algorithm import Algorithm


class GCSQ(Algorithm):
    def __init__(self, seed):
        super().__init__(seed)
        self.name = "GCS-Q"

    # implementation is a tiny bit different from the one in paper, but does the same thing
    def solve(self, n, edges, num_coalitions=None):
        # initialize with grand coalition
        coalitions = [list(range(n))]
        # for a maximum of n steps (as we can at most have n coalitions, one for each agent)
        for _ in range(n):
            new_coalitions = copy.deepcopy(coalitions)
            for c in coalitions:
                # unless the coalition cannot be further split because it already contains only one agent
                if len(c) > 1:
                    # Solve optimal split problem
                    # Let annealer find optimal split of coalition
                    c1, c2 = self._split(c, edges)
                    # If splitting increases value of this part of the coalition structure,
                    # replace coalition c with split coalitions
                    if utils.value(c1, edges) + utils.value(c2, edges) > utils.value(c, edges):
                        new_coalitions.remove(c)
                        new_coalitions.append(c1)
                        new_coalitions.append(c2)
            # none of the splittings of the coalitions of the previous coalition structure improved it, stop the algorithm
            # and return the previous coalition structure
            if len(coalitions) == len(new_coalitions):
                break
            coalitions = new_coalitions
        return coalitions

    def _split(self, c, edges):
        Q = {}
        for i in range(len(c)):
            for j in range(len(c)):
                if i < j:
                    utils.add(Q, i, i, edges[(c[i],c[j])])
                    utils.add(Q, j, j, edges[(c[i],c[j])])
                    utils.add(Q, i, j, -2*edges[(c[i],c[j])])
        solution = utils.solve_with_qbsolv(Q, len(c), self.seed)
        c1 = [c[k] for k in range(len(c)) if solution[k] == 1]
        c2 = [c[k] for k in range(len(c)) if solution[k] == 0]
        return c1, c2


