# Code Grade Challenge, UCSD Genome Course Week2
# Find Eulerian Cycle from adjacency list
# Input: Adjacency list kmers separated by "->"
# Return: Eulerian path through all kmers
# Created July 2018, Callen Hyland


import sys
from collections import defaultdict


def EulerianPathFromAdjacency(ad):
    nodes = ad.keys()
    paths = dict.fromkeys(nodes, 0)
    for node in nodes:
        paths[node] += len(ad[node])
        for out in ad[node]:
            if out in paths:
                paths[out] -= 1
            else:
                paths[out] = -1

    current_location = [i for i in paths if  paths[i] == 1][0]
    stack = []
    circuit = []

    while bool(ad):
        while ad[current_location] != []:
            stack.append(current_location)
            current_location = ad[current_location].pop()
        circuit.append(current_location)
        if stack != []:
            current_location = stack.pop()
        else:
            return(circuit[::-1])


adj_dict = defaultdict(list)
lines = sys.stdin.read().splitlines()
for line in lines:
    paths = line.rstrip().split(" -> ")
    adj_dict[paths[0]] = paths[1].split(",")

circuit = EulerianPathFromAdjacency(adj_dict)
print("->".join(circuit))

