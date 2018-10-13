# -*- coding: utf-8 -*-
"""
Code Challenge: Solve the String Reconstruction Problem.

 Input: An integer k followed by a list of k-mers Patterns.
 Output: A string Text with k-mer composition equal to Patterns. 
 (If multiple answers exist, you may return any one.)
 
 Created by Callen Hyland, August 2018
"""

from collections import defaultdict


def EulerianPathFromAdjacency(ad):
    
    # Find location with the extra outgoing node
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
    
    #initialize stack and circuit
    stack = []
    circuit = []
    
    # Recreate Eulerian circuit through graph
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
    # Create a DeBruijn graph from a list of k-mers
    db = defaultdict(list)
    for kmer in kmers:
        db[kmer[0:-1]].append(kmer[1:])
    return(db)


def assemble(lines):
    # Assemble genome from list of overlapping lines
    genome = lines[0]
    for i in range(1,len(lines)):
        genome = genome + lines[i][-1]
    return(genome)


"""
Script for reading file, creating Eulerian path from adjacency list,
and assembling into a single string ("genome")
""" 

# Read lines from text file as kmer list
lines = []
with open ("dataset_203_7.txt", "r") as myfile:
    for line in myfile:
        lines.append(line.rstrip())
kmers = lines[1:]

db = debruijn_from_kmers(kmers)
path = EulerianPathFromAdjacency(db)
genome = assemble(path)

# Write genome sequence to output file
f = open("genome-from-kmers-step7.txt", "w")
f.write(genome)
f.close()
