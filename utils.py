from dwave_qbsolv import QBSolv
import numpy as np



# this function solves a given QUBO-Matrix Q with Qbsolv
def solve_with_qbsolv(Q, n, seed):
    response = QBSolv().sample_qubo(Q, num_repeats=100, timeout=5, seed=seed)
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



# this function prints the first n row/columns of a QUBO-Matrix Q
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



def generate_problem(n, mean=0.6):
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
