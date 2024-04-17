import networkx as nx
import pickle
import numpy as np

def transform_for_nx(naeimeh_graph):
    edges = []
    for key in naeimeh_graph.keys():
        i, j = key.split(",")
        i, j = int(i)-1, int(j)-1
        edges.append((i, j, naeimeh_graph[key]))
    return edges

data = pickle.load(open('data_new_20samples_4_28.pkl', 'rb'))

graph_sizes = list(data.keys())  # is [4,6,8,10,12,14,16,18,20,22,24,26,28] for E.ON data
num_graphs_per_size = len(data[graph_sizes[0]])  # is 20 for E.ON data
num_graph_sizes = len(graph_sizes)
for num_agents in graph_sizes:
    for graph_num in range(num_graphs_per_size):
        graph = data[num_agents][graph_num]
        edges = transform_for_nx(graph)
        G = nx.Graph()
        G.add_weighted_edges_from(edges)
        nx.write_weighted_edgelist(G, f"eon_data_edgelists/eon_graph_size_{num_agents}_num_{graph_num}.edgelist")

# Also, while we are at it anyway, determine some stats on the graph weights
max = 0
min = 0
zero_weighted = 0
positive = 0
negative = 0
list = []
for key, graph_size_list in data.items():
    for graph_list in graph_size_list:
        for key, value in graph_list.items():
            if value < min:
                min = value
            elif value > max:
                max = value
            list.append(value)
            if value > 0:
                positive = positive + 1
            elif value < 0:
                negative = negative + 1
            elif value == 0:
                zero_weighted = zero_weighted + 1
print(f"Max: {max}, Min: {min}")
list = np.array(list)
avg = np.mean(list)
std = np.std(list)
print(f"Average: {avg}, Std: {std}")
print(f"Zero-weighted: {zero_weighted}, positive: {positive}, negative: {negative}; non-zero-weighted: {positive + negative}")