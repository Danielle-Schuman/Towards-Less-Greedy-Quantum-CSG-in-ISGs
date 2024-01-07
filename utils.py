import copy
from dwave_qbsolv import QBSolv
import numpy as np



# this function solves a given QUBO-Matrix Q with Qbsolv
def solve_with_qbsolv(Q, n, seed, timeout=10):
    response = QBSolv().sample_qubo(Q, num_repeats=1000, timeout=timeout, seed=seed)
    solution = [response.samples()[0][i] for i in range(n)]
    return solution



# this function calculates the value of a solution for a given QUBO-Matrix Q
def getValue(Q, solution):
    ones = [i for i in range(len(solution)) if solution[i] == 1]
    value = 0
    for x in ones:
        for y in ones:
            if (x,y) in Q.keys():
                value += Q[(x,y)]
    return value


def transform(naeimeh_graph):
    edges = {}
    for key in naeimeh_graph.keys():
        i, j = key.split(",")
        i, j = int(i)-1, int(j)-1
        edges[(i,j)] = naeimeh_graph[key]
    return edges


# this function prints the first num_agents row/columns of a QUBO-Matrix Q
def printQUBO(Q, n):
    for row in range(n):
        for column in range(n):
            if row > column:
                print("        ", end = '')
                continue
            printing = ""
            if (row,column) in Q.keys() and Q[(row,column)] != 0:
                printing = str(float(Q[(row,column)]))
            printing += "_______"
            printing = printing[:7]
            printing += " "
            print(printing, end = '')
        print("")



def generate_problem(n, mean=0.5):
    edges = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                edges[(i,j)] = np.random.sample() - mean
    return edges



def add(Q, i, j, v):
    if (i,j) not in Q.keys():
        Q[(i,j)] = v
    else:
        Q[(i,j)] += v



def value(c, edges):
    v = 0
    for i in c:
        for j in c:
            if i < j:
                v += edges[(i,j)]
    return v


def baseline(n, edges):

    best_coalitions = None
    best_value = 0

    for _ in range(10):
        agent_list = np.array(range(n))
        np.random.shuffle(agent_list)
        coalitions = sub_baseline(edges, agent_list)
        v = np.sum([value(c,edges) for c in coalitions])
        if v > best_value:
            best_coalitions = coalitions
            best_value = v

    return best_coalitions


def sub_baseline(edges, agent_list):
    coalitions = [[agent_list[0]]]
    for i in agent_list[1:]:
        improvements = []
        for c in coalitions:
            old_v = value(c, edges)
            new_c = copy.deepcopy(c)
            new_c.append(i)
            new_v = value(new_c, edges)
            improvements.append(new_v - old_v)
        if np.max(improvements) > 0:
            coalitions[np.argmax(improvements)].append(i)
        else:
            coalitions.append([i])
    return coalitions