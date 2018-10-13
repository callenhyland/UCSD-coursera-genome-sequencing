# -*- coding: utf-8 -*-
"""
Functions for Week 2 of UCSD Genome Sequencing Course
Callen Hyland, August 2018
"""

from collections import defaultdict

def EulerianCycleFromAdjacency(adj):
    #choose current location and initialize stack and circuit
    stack = []
    circuit = []
    current_location = "0"

    while bool(adj):
        while adj[current_location] != []:
            stack.append(current_location)
            current_location = adj[current_location].pop()
        circuit.append(current_location)
        if stack != []:
            current_location = stack.pop()
        else:
            return(circuit[::-1])


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


def debruijn_from_kmers(kmers):
    db = defaultdict(list)
    for kmer in kmers:
        db[kmer[0:-1]].append(kmer[1:])
    #nodes = sorted(db.keys())
    #return([(node + " -> " + ",".join(db[node])) for node in nodes])
    return(db)
    

def assemble(lines):
    genome = lines[0]
    for i in range(1,len(lines)):
        genome = genome + lines[i][-1]
    return(genome)






