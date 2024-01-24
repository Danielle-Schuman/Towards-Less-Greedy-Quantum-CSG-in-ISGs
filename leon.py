import utils

from algorithm import QuantumAlgorithm

def _find_region_edges(region_list):
    region_edges_list =[]
    for region in region_list:
        region_edges = []
        for i in range(len(region)-1):
            region_edges.append((region[i], region[i+1]))
        region_edges.append((region[-1], region[0]))
        region_edges_list.append(region_edges)
    return region_edges_list

def _star_analysis(feeder_graph, region, region_list, edges_to_add):
    region_list_edges = _find_region_edges(region_list)
def _planarize_old(num_agents, input_graph):
    # Step 0: Create minimal input graph K4
    feeder_graph = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
    region_list = [[0, 1, 3], [1, 2, 3], [2, 0, 3], [0, 2, 1]]
    vertex_to_add = 4
    while vertex_to_add < num_agents:
        # Step 1: Construct list of edges to add to feeder graph from vertex_to_add to all other vertices
        edges_to_add = [(i, vertex_to_add) for i in range(vertex_to_add)]
        # Step 2: Use Star Analysis to generate Path List descriptions for minimal Kn for each region
        for region in range(len(region_list)):  # TODO parallelize this
            _star_analysis(feeder_graph, region, region_list, edges_to_add)
        # Step 3: Select globally minimal graphs Kn
        # Step 4: Create graph descriptions for all selected graphs
        # Step 5: Format conversion for isomorphic testing
        # Step 6: Categorization into isomorphic families by pair-wise isomorphic testing
        # Step 7: Chose one representative from each isomorphic family
        pass  # TODO remove this line

def _planarize(num_agents, input_graph):
    pass

class Leon(QuantumAlgorithm):
    def __init__(self, seed, num_graph_sizes, solver="qbsolv", num_coalitions=None, timeout=600):
        super().__init__(seed=seed, num_graph_sizes=num_graph_sizes, num_coalitions=num_coalitions, solver=solver, timeout=timeout)
        self.name = f"Leon19_{self.solver}"

    def solve(self, num_agents, edges):
        #TODO
        # Step 1: Planarize
        edges = _planarize(num_agents, edges)
        # Step 2: Create logical qubit lattice (see Fig. 7)
        # Step 3: Create physical qubit lattice
        # Step 3.5: Convert physical qubit lattice to qubo
        pass
        