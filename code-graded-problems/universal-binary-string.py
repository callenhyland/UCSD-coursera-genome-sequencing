# -*- coding: utf-8 -*-
"""
Code Challenge: N-univeral circular binary string problem

 Input: An integer n, the length of binary substrings
 Output: reconstructed N-universal circular binary string
 
 Created by Callen Hyland, August 2018
"""

from collections import defaultdict
from itertools import product
from random import choice


def debruijn_from_kmers(kmers):
    db = defaultdict(list)
    for kmer in kmers:
        db[kmer[0:-1]].append(kmer[1:])
    return(db)


def EulerianCycleFromAdjacency(adj):
    #choose current location and initialize stack and circuit
    stack = []
    circuit = []
    current_location = choice(list(adj.keys()))

    while bool(adj):
        while adj[current_location] != []:
            stack.append(current_location)
            current_location = adj[current_location].pop()
        circuit.append(current_location)
        if stack != []:
            current_location = stack.pop()
        else:
            return(circuit[::-1])

        
def assemble(lines):
    genome = lines[0]
    for i in range(1,len(lines)):
        genome = genome + lines[i][-1]
    return(genome)


"""
Solve the N-Universal Circular Binary String problem
Inputs: n (length of binary sub-strings)
Output: prints reconstructed string
"""

n = 8 # length of binary strings
bn = ["".join(seq) for seq in product("01", repeat = n)]

db = debruijn_from_kmers(bn)
cycle = EulerianCycleFromAdjacency(db)
print(assemble(cycle[:-(n-1)]))